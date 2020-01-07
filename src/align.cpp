/*
 * Some code used for upgrading Sapling to a full read aligner
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

static void write_sam_alignment(StripedSmithWaterman::Alignment* a, string read_name, string read_qual, string read_seq, string ref_name, int strand, int aligned, FILE *fout)
{
  fprintf(fout, "%s\t", read_name.c_str());
	if (!aligned) fprintf(fout, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
	else
  {
    int32_t c, p;
    uint32_t mapq = -4.343 * log(1 - (double)abs(a->sw_score - a->sw_score_next_best)/(double)a->sw_score);
    mapq = (uint32_t) (mapq + 4.99);
    mapq = mapq < 254 ? mapq : 254;
    if (strand) fprintf(fout, "16\t");
    else fprintf(fout, "0\t");
    fprintf(fout, "%s\t%d\t%d\t", ref_name.c_str(), a->ref_begin + 1, mapq);
    for (c = 0; c < (int)(a->cigar.size()); ++c)
    {
	    char letter = cigar_int_to_op(a->cigar[c]);
	    uint32_t length = cigar_int_to_len(a->cigar[c]);
	    fprintf(fout, "%lu%c", (unsigned long)length, letter);
    }
		fprintf(fout, "\t*\t0\t0\t");
		fprintf(fout, "%s", read_seq.c_str());
		fprintf(fout, "\t");

    if (read_qual.length() > 0 && strand)
    {
	    for (p = read_qual.length() - 1; p >= 0; --p) fprintf(fout, "%c", read_qual[p]);
    }
    else if (read_qual.length() > 0)
    {
      fprintf (fout, "%s", read_qual.c_str());
    }
    else
    {
      fprintf(fout, "*");
    }
    fprintf(fout, "\tAS:i:%d", a->sw_score);
    fprintf(fout,"\tNM:i:%d\t", a->mismatches);
    if (a->sw_score_next_best > 0) fprintf(fout, "ZS:i:%d\n", a->sw_score_next_best);
    else fprintf(fout, "\n");
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
            string refString = string(refFn);
            sapling = new Sapling(refString, refString + "_k16.sa", refString + "_k16.sap", -1, -1, 16, "");
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
        void align_all_reads(char *outFn)
        {
            cout << "Aligning reads" << endl;
            FILE *fout = fopen(outFn, "w");
            fprintf(fout, "@HD\tVN:1.6\tSO:coordinate\n");
            size_t lastEnd = 0;
            for (map<string, int>::iterator it = sapling->chrEnds.begin(); it != sapling->chrEnds.end(); it++ )
            {
              
              size_t endCoord = it->second;
              size_t chrLength = endCoord - lastEnd;
              string name = it->first;
              fprintf(fout, "@SQ\tSN:%s\tLN:%zu\n", name.c_str(), chrLength);
              lastEnd = endCoord;
            }
            while(1) 
            {
                string cur = align_next_read(fout);
                if(cur.length() == 0)
                {
                    return;
                }
                //cout << cur << endl;
            }
            fclose(fout);
        }

        // Align the next read in the stream and return an alignment record
        string align_next_read(FILE *fout)
        {
            vector<string> cur_read = get_read();
            if(cur_read.size() < 3)
            {
                return "";
            }
            string seq = cur_read[1];
            string name = cur_read[0].substr(1);
            string qual = cur_read[2];
            return seed_extend(seq, name, qual, fout);
        }

        char complement(char c)
        {
            if(c == 'A') return 'T';
            if(c == 'C') return 'G';
            if(c == 'G') return 'C';
            if(c == 'T') return 'A';
            return c;
        } 

        string revComp(string s)
        {
          string res(s.length(), 'A');
          for(size_t i = 0; i<s.length(); i++) res[i] = complement(s[s.length()-1-i]);
          return res;
        }

        string seed_extend(string &readSeq, string &name, string &qual, FILE *fout)
        {
            size_t num_seeds = 7;
            size_t last = readSeq.length() - sapling->k;
            int best_score = -1;
            int bestStrand = 0;
            StripedSmithWaterman::Alignment best_alignment;
            for(int iter = 0; iter < 2; iter++)
            {
              string seq = iter ? revComp(readSeq) : readSeq;
              vector<tuple<size_t, size_t, size_t, size_t, size_t>> counts;
              //cout << "Iter: " << iter << " " << seq << endl;
              for(size_t i = 0; i < num_seeds; i++)
              {
                  size_t cur_pos = 0;
                  if(i == num_seeds - 1) cur_pos = last;
                  else if(i > 0) cur_pos = last / (num_seeds - 1) * i;

                  string query = seq.substr(cur_pos, sapling->k);
                  long long val = sapling->kmerize(query);
                  long long ref_pos_signed = sapling->plQuery(query, val, sapling->k);
                  size_t ref_pos = 0;
                  if(ref_pos_signed == -1) continue;
                  else ref_pos = (size_t)ref_pos_signed;
                  string ref_seq = sapling->reference.substr(ref_pos, sapling->k);
                  int strcmp = query.compare(ref_seq);
                  //cout << query << " " << ref_seq << endl;  
                  if(!strcmp)
                  {
                      //cerr << "Getting range for match at position " << ref_pos << endl;
                      //cerr << "Sequence = " << query << endl;
                      size_t sa_pos = sapling->sa[ref_pos];
                      int left = sapling->countHitsLeft(sapling->sa[ref_pos]);
                      int right = sapling->countHitsRight(sapling->sa[ref_pos]);
                      
                      //left = min(left, (int)(sapling->MAX_HITS/2));
                      //right = min(right, (int)(sapling->MAX_HITS/2));
                      counts.push_back(
                        make_tuple(
                          (size_t)(left + right + 1),
                          cur_pos,
                          sa_pos,
                          left,
                          right)
                        );
                  }
                  else
                  {
                      //cout << "No matches for seed " << i << " of " << name << " found" << endl;
                  }
              }
              sort(counts.begin(), counts.end());
              for(size_t i = 0; i<counts.size(); i++)
              {
                  size_t query_pos = get<1>(counts[i]);
                  size_t left_flank = query_pos;
                  size_t right_flank = seq.length() - query_pos;
                  size_t sa_pos = get<2>(counts[i]);
                  int left = get<3>(counts[i]); 
                  int right = get<4>(counts[i]);
                  if((size_t)(left + right) > sapling->MAX_HITS)
                  {
                    if(best_score == -1)
                    {
                      left = min(left, (int)(sapling->MAX_HITS/2));
                      right = min(right, (int)(sapling->MAX_HITS/2));
                    }
                    else
                    {
                      left = right = 0;
                    }
                  }

                  for(int offset = -left; offset <= right; offset++)
                  {
                      size_t ref_pos = sapling->rev[sa_pos + offset];
                      int ref_left_pos = (int)ref_pos - left_flank - 10;
                      if(ref_left_pos < 0) ref_left_pos = 0;
                      size_t ref_right_pos = (int)ref_pos + right_flank + 10;
                      if(ref_right_pos >= sapling->n) ref_right_pos = sapling->n;
                      int ref_length = ref_right_pos - ref_left_pos;
                      string ref_seq = sapling->reference.substr(ref_left_pos, ref_length);
                      StripedSmithWaterman::Alignment aln_result;
                      StripedSmithWaterman::Filter filter;
                      aln->Align(seq.c_str(), ref_seq.c_str(),
                              ref_seq.length(), filter, &aln_result, 15);
                      uint16_t cur_score = aln_result.sw_score;
                      aln_result.ref_begin += ref_left_pos;
                      if((int)cur_score > best_score)
                      {
                          best_score = cur_score;
                          best_alignment = aln_result;
                          bestStrand = iter;
                      }
                  }    
              }
            }
            //cout << best_score << " " << bestStrand << endl;
            if(best_score > -1)
            {
                string refName = "*";
                size_t bestEnd = 0;
                size_t lastEnd = 0;
                map<string, int>::iterator it;

                for ( it = sapling->chrEnds.begin(); it != sapling->chrEnds.end(); it++ )
                {
                  size_t endCoord = it->second;
                  if((long long)endCoord > best_alignment.ref_begin && (bestEnd == 0 || endCoord < bestEnd))
                  {
                    bestEnd = endCoord;
                    refName = it->first;
                  }
                  if((long long)endCoord <= best_alignment.ref_begin && (lastEnd == 0 || endCoord > lastEnd))
                  {
                    lastEnd = endCoord;
                  }
                }
                best_alignment.ref_begin -= lastEnd;
                write_sam_alignment(&best_alignment, name, qual, readSeq, refName, bestStrand, 1, fout);
            }
            else
            {
              write_sam_alignment(&best_alignment, name, qual, readSeq, "refname", bestStrand, 0, fout);
            }
            return "success";
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
    printf("usage: ./align <query> <ref> <outfile>\n");
}

int main(int argc, char** argv)
{
    if(argc != 4)
    {
        usage();
        return 1;
    }
    char* queryFn = argv[1];
    char* refFn = argv[2];
    char* outFn = argv[3];
    SaplingAligner* sa = new SaplingAligner(queryFn, refFn);
    sa->align_all_reads(outFn);
    delete sa;
    return 0;
}
