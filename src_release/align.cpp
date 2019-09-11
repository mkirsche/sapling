/*
 * A work in progress upgrading Sapling to a full read aligner
 */
#include "ssw_cpp.h"
#include "sapling_api.h"
#include <fstream>
#include <sstream>
#include <string>
#include <bits/stdc++.h> 

#define MAPSTR "MIDNSHP=X"
#ifndef BAM_CIGAR_SHIFT
#define BAM_CIGAR_SHIFT 4u
#endif

static inline char cigar_int_to_op(uint32_t cigar_int) {
	return (cigar_int & 0xfU) > 8 ? 'M': MAPSTR[cigar_int & 0xfU];
}

/*!	@function		Extract length of a CIGAR operation from CIGAR 32-bit unsigned integer
	@param	cigar_int	32-bit unsigned integer, representing encoded CIGAR operation and length
	@return			length of CIGAR operation
*/
static inline uint32_t cigar_int_to_len (uint32_t cigar_int) {
	return cigar_int >> BAM_CIGAR_SHIFT;
}

// Printing an alignment for debugging purposes
static void PrintAlignment(const StripedSmithWaterman::Alignment& alignment){
  cout << "===== SSW result =====" << endl;
  cout << "Best Smith-Waterman score:\t" << alignment.sw_score << endl
       << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << endl
       << "Reference start:\t" << alignment.ref_begin << endl
       << "Reference end:\t" << alignment.ref_end << endl
       << "Query start:\t" << alignment.query_begin << endl
       << "Query end:\t" << alignment.query_end << endl
       << "Next-best reference end:\t" << alignment.ref_end_next_best << endl
       << "Number of mismatches:\t" << alignment.mismatches << endl
       << "Cigar: " << alignment.cigar_string << endl;
  cout << "======================" << endl;
}

static void write_sam_alignment(StripedSmithWaterman::Alignment* a, string read_name, string read_qual, string read_seq, string ref_name, int strand)
{
    fprintf(stdout, "%s\t", read_name.c_str());
	if (a->sw_score == 0) fprintf(stdout, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
	else {
		int32_t c, p;
		uint32_t mapq = -4.343 * log(1 - (double)abs(a->sw_score - a->sw_score_next_best)/(double)a->sw_score);
		mapq = (uint32_t) (mapq + 4.99);
		mapq = mapq < 254 ? mapq : 254;
		if (strand) fprintf(stdout, "16\t");
		else fprintf(stdout, "0\t");
		fprintf(stdout, "%s\t%d\t%d\t", ref_name.c_str(), a->ref_begin + 1, mapq);
		for (c = 0; c < (int)(a->cigar.size()); ++c)
        {
			char letter = cigar_int_to_op(a->cigar[c]);
			uint32_t length = cigar_int_to_len(a->cigar[c]);
			fprintf(stdout, "%lu%c", (unsigned long)length, letter);
		}
		fprintf(stdout, "\t*\t0\t0\t");
		fprintf(stdout, "%s", read_seq.c_str());
		fprintf(stdout, "\t");
		if (read_qual.length() > 0 && strand) {
			for (p = read_qual.length() - 1; p >= 0; --p) fprintf(stdout, "%c", read_qual[p]);
		}else if (read_qual.length() > 0) fprintf (stdout, "%s", read_qual.c_str());
		else fprintf(stdout, "*");
		fprintf(stdout, "\tAS:i:%d", a->sw_score);
		fprintf(stdout,"\tNM:i:%d\t", a->mismatches);
		if (a->sw_score_next_best > 0) fprintf(stdout, "ZS:i:%d\n", a->sw_score_next_best);
		else fprintf(stdout, "\n");
    }
}

/*
 * The structure containing the high-level alignment logic
 */
class SaplingAligner
{
    // A dynamic programming aligner for the extension
    StripedSmithWaterman::Aligner* aln;

    // A binary search-based seed finder
    Sapling* sapling;

    // A stream which inputs reads to process
    std::ifstream *queryReader;

    public:

        // Basic constructor which initializes everything
        SaplingAligner(char* queryFn, char* refFn)
        {
            sapling = new Sapling(refFn);
            aln = new StripedSmithWaterman::Aligner();
            queryReader = new std::ifstream(queryFn);
        }

        // Get the next read from the input stream
        // For now just return the name and sequence
        vector<string> get_read()
        {
            vector<string> res;
            string line;
            for(size_t i = 0; i<4; i++)
            {
                if(!std::getline(*queryReader, line))
                {
                    return res;
                }
                if(i == 0 || i == 1 ||i == 3)
                {
                    res.push_back(line);
                }
            }
            return res;
        }

        // Run the alignment of all reads in the input file
        void align_all_reads()
        {
            while(1) 
            {
                string cur = align_next_read();
                if(cur.length() == 0)
                {
                    return;
                }
                cout << cur << endl;
            }
        }

        // Align the next read in the stream and return an alignment record
        string align_next_read()
        {
            vector<string> cur_read = get_read();
            if(cur_read.size() < 3)
            {
                return "";
            }
            string seq = cur_read[1];
            string name = cur_read[0].substr(1);
            string qual = cur_read[2];
            return seed_extend(seq, name, qual);
        }

        string seed_extend(string seq, string name, string qual)
        {
            size_t num_seeds = 3;
            size_t last = seq.length() - sapling->k;
            vector<tuple<size_t, size_t, size_t, size_t, size_t>> counts;
            for(size_t i = 0; i < num_seeds; i++)
            {
                size_t cur_pos = 0;
                if(i == num_seeds - 1) cur_pos = last;
                else if(i > 0) cur_pos = last / (num_seeds - 1) * i;

                string query = seq.substr(cur_pos, sapling->k);
                long long val = sapling->kmerize(query);
                size_t ref_pos = sapling->plQuery(query, val, sapling->k);
                string ref_seq = sapling->reference.substr(ref_pos, sapling->k);
                int strcmp = query.compare(ref_seq);
                if(!strcmp)
                {
                    cout << "Getting range for match at position " << ref_pos << endl;
                    cout << "Sequence = " << query << endl;
                    size_t sa_pos = sapling->sa[ref_pos];
                    int left = sapling->countHitsLeft(sapling->sa[ref_pos]);
                    int right = sapling->countHitsRight(sapling->sa[ref_pos]);
                    if((size_t)(left + right) < sapling->MAX_HITS)
                    {
                        counts.push_back(
                            make_tuple(
                                (size_t)(left + right + 1),
                                cur_pos,
                                sa_pos,
                                left,
                                right)
                        );
                        for(int offset = -left; offset<=right; offset++)
                        {
                            //size_t real_pos = sapling->rev[sa_pos + offset];
                            //cout << "Match of seed " << i << " at " << real_pos << endl;
                            //cout << sapling->reference.substr(real_pos, sapling->k) << endl;
                        }
                    }
                }
                else
                {
                    //cout << "No matches for seed " << i << " found" << endl;
                }
            }
            sort(counts.begin(), counts.end());
            int best_score = -1;
            StripedSmithWaterman::Alignment best_alignment;
            for(size_t i = 0; i<counts.size(); i++)
            {
                size_t query_pos = get<1>(counts[i]);
                size_t left_flank = query_pos;
                size_t right_flank = seq.length() - query_pos;
                size_t sa_pos = get<2>(counts[i]);
                int left = get<3>(counts[i]);
                int right = get<4>(counts[i]);


                for(int offset = -left; offset <= right; offset++)
                {
                    size_t ref_pos = sapling->rev[sa_pos + offset];
                    //cout << query_pos << " " << ref_pos << endl;
                    int ref_left_pos = (int)ref_pos - left_flank - 10;
                    if(ref_left_pos < 0) ref_left_pos = 0;
                    size_t ref_right_pos = (int)ref_pos + right_flank + 10;
                    if(ref_right_pos >= sapling->n) ref_right_pos = sapling->n;
                    int ref_length = ref_right_pos - ref_left_pos;
                    string ref_seq = sapling->reference.substr(ref_left_pos, ref_length);
                    //cout << seq << " " << ref_seq << endl; 
                    StripedSmithWaterman::Alignment aln_result;
                    StripedSmithWaterman::Filter filter;
                    aln->Align(seq.c_str(), ref_seq.c_str(),
                            ref_seq.length(), filter, &aln_result, 15);
                    uint16_t cur_score = aln_result.sw_score;
                    if((int)cur_score > best_score)
                    {
                        best_score = cur_score;
                        best_alignment = aln_result;
                    }
                }    
            }
            if(best_score > -1)
            {
                //PrintAlignment(best_alignment);
                write_sam_alignment(&best_alignment, name, qual, seq, "refname", 0);
            }
            return "";
        }

        // Deconstructor
        ~SaplingAligner()
        {
            delete aln;
            delete sapling;
            delete queryReader;
        }
};

// Error message for when the user provides the wrong parameters
void usage()
{
    printf("usage: ./align <query> <ref>\n");
}

int main(int argc, char** argv)
{
    if(argc != 3)
    {
        usage();
        return 1;
    }
    char* queryFn = argv[1];
    char* refFn = argv[2];
    SaplingAligner* sa = new SaplingAligner(queryFn, refFn);
    sa->align_all_reads();
    delete sa;
    return 0;
}
