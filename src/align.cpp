/*
 * A seed-and-extend aligner using Sapling to augment a suffix array for fast seeding
 * and a striped-smith-waterman algorithm for performing extension
 */
#include "ssw_cpp.h"
#include "sapling_api.h"
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#define MAPSTR "MIDNSHP=X"
#ifndef BAM_CIGAR_SHIFT
#define BAM_CIGAR_SHIFT 4u
#endif

char* queryFn;
char* refFn;
char* outFn;
size_t num_seeds = 7;
size_t flankingSequence = 2;
int saplingK = 16;
size_t maxHits = 32;

/*
 * Error message for when the user provides the wrong parameters
 */
void usage()
{
    printf("usage: ./align <query> <ref> <outfile> [num_seeds=<int>] [sapling_k=<int>] [flanking_sequence=<int>] [max_hits=<int>]\n");
}

/*
 * Parser for the command line arguments
 */
void parseArgs(int argc, char** argv)
{
  queryFn = argv[1];
  refFn = argv[2];
  outFn = argv[3];
  for(int i = 4; i<argc; i++)
  {
    string cur = argv[i];
    size_t eqPos = cur.find("=");
    if(eqPos != string::npos)
    {
      string arg = cur.substr(0, eqPos);
      string val = cur.substr(eqPos + 1);
      if(arg.compare("num_seeds") == 0)
      {
        num_seeds = stoi(val);
      }
      else if(arg.compare("sapling_k") == 0)
      {
        saplingK = stoi(val);
      }
      else if(arg.compare("flanking_sequence") == 0)
      {
        flankingSequence = stoi(val);
      } 
      else if(arg.compare("max_hits") == 0)
      {
        maxHits = stoi(val);
      }  
    }
  }
}

/*
 * Function for mapping integers representing CIGAR operations to their single-character abbreviations
 */
static inline char cigar_int_to_op(uint32_t cigar_int) {
	return (cigar_int & 0xfU) > 8 ? 'M': MAPSTR[cigar_int & 0xfU];
}

/*
 * Extract length of a CIGAR operation from CIGAR 32-bit unsigned integer
 */
static inline uint32_t cigar_int_to_len (uint32_t cigar_int) {
	return cigar_int >> BAM_CIGAR_SHIFT;
}

/*
 * Writes an alignment to a file in SAM format
 */
static void write_sam_alignment(StripedSmithWaterman::Alignment* a, string read_name, string read_qual, string read_seq, string ref_name, int strand, int aligned, FILE *fout)
{
  fprintf(fout, "%s\t", read_name.c_str());

  // Unaligned read case
	if (!aligned)
  {
    fprintf(fout, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
  }

  // Valid alignment so output fields one-by-one
	else
  {
    int32_t c, p;

    // Compute mapping quality
    uint32_t mapq = -4.343 * log(1 - (double)abs(a->sw_score - a->sw_score_next_best)/(double)a->sw_score);
    mapq = (uint32_t) (mapq + 4.99);
    mapq = mapq < 254 ? mapq : 254;

    // Output 0 if forward strand and 16 if reverse
    if (strand) fprintf(fout, "16\t");
    else fprintf(fout, "0\t");

    // Print the reference name, position, and mapping quality
    fprintf(fout, "%s\t%d\t%d\t", ref_name.c_str(), a->ref_begin + 1, mapq);

    // Print the CIGAR string for the alignment
    for (c = 0; c < (int)(a->cigar.size()); ++c)
    {
	    char letter = cigar_int_to_op(a->cigar[c]);
	    uint32_t length = cigar_int_to_len(a->cigar[c]);
	    fprintf(fout, "%lu%c", (unsigned long)length, letter);
    }

    // Print the read sequence
		fprintf(fout, "\t*\t0\t0\t");
		fprintf(fout, "%s", read_seq.c_str());
		fprintf(fout, "\t");

    // Print the read quality information
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

    // Print the alignment score, mismatch count, and next best score
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
            sapling = new Sapling(refString, refString + ".sa", refString + "_k" + to_string(saplingK) + ".sap", -1, -1, saplingK, "");
            aln = new StripedSmithWaterman::Aligner();
            queryReader = new std::ifstream(queryFn);
        }

        // Get the next read from the input stream
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
        void align_all_reads(char *outFn, int argc, char** argv)
        {
            cout << "Aligning reads" << endl;
            FILE *fout = fopen(outFn, "w");
            fprintf(fout, "@HD\tVN:1.6\tSO:coordinate\n");
            size_t lastEnd = 0;
            for (map<size_t, string>::iterator it = sapling->chrEnds.begin(); it != sapling->chrEnds.end(); it++ )
            {
              
              size_t endCoord = it->first;
              size_t chrLength = endCoord - lastEnd;
              string name = it->second;
              fprintf(fout, "@SQ\tSN:%s\tLN:%zu\n", name.c_str(), chrLength);
              lastEnd = endCoord;
            }
            fprintf(fout, "@PG\tID:sapling\tVN:1.0\tCL:%s", argv[0]);
            for(int i = 1; i<argc; i++)
            {
              fprintf(fout, " %s", argv[i]);
            }
            fprintf(fout, "\n");
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

        // The complement of a given base-pair character
        char complement(char c)
        {
            if(c == 'A') return 'T';
            if(c == 'C') return 'G';
            if(c == 'G') return 'C';
            if(c == 'T') return 'A';
            return c;
        } 

        // The reverse complement of a genomic string
        string revComp(string s)
        {
          string res(s.length(), 'A');
          for(size_t i = 0; i<s.length(); i++) res[i] = complement(s[s.length()-1-i]);
          return res;
        }

        // Performs seed-and-extend alignment of a given read and writes the alignment to a file
        string seed_extend(string &readSeq, string &name, string &qual, FILE *fout)
        {
            size_t last = readSeq.length() - sapling->k;
            int best_score = -1;
            int bestStrand = 0;
            StripedSmithWaterman::Alignment best_alignment;
            size_t best_offset = 0;
            bool done = 0;
            for(int iter = 0; iter < 2 && !done; iter++)
            {
              string seq = iter ? revComp(readSeq) : readSeq;
              vector<tuple<size_t, size_t, size_t, size_t, size_t>> counts;
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
                  if(!strcmp)
                  {
                      size_t sa_pos = sapling->sa[ref_pos];
                      size_t left = sapling->countHitsLeft(sapling->sa[ref_pos], maxHits);
                      size_t right = sapling->countHitsRight(sapling->sa[ref_pos], maxHits);

                      counts.push_back(
                        make_tuple(
                          (left + right + 1),
                          cur_pos,
                          sa_pos,
                          left,
                          right)
                        );
                  }
              }
              sort(counts.begin(), counts.end());
              for(size_t i = 0; i<counts.size() && !done; i++)
              {
                  size_t query_pos = get<1>(counts[i]);
                  size_t left_flank = query_pos;
                  size_t right_flank = seq.length() - query_pos;
                  size_t sa_pos = get<2>(counts[i]);
                  int left = get<3>(counts[i]); 
                  int right = get<4>(counts[i]);
                  if((size_t)(left + right) > maxHits)
                  {
                    if(best_score == -1)
                    {
                      left = min(left, (int)(maxHits/2));
                      right = min(right, (int)(maxHits/2));
                    }
                    else
                    {
                      left = right = 0;
                    }
                  }

                  for(int offset = -left; offset <= right && !done; offset++)
                  {
                      size_t ref_pos = sapling->rev[sa_pos + offset];
                      long long ref_left_pos = (long long)ref_pos - left_flank - flankingSequence;
                      if(ref_left_pos < 0) ref_left_pos = 0;
                      long long ref_right_pos = (long long)ref_pos + right_flank + flankingSequence;
                      if((size_t)ref_right_pos >= sapling->n) continue;
                      int ref_length = ref_right_pos - ref_left_pos;
                      string ref_seq = sapling->reference.substr(ref_left_pos, ref_length);
                      StripedSmithWaterman::Alignment aln_result;
                      StripedSmithWaterman::Filter filter;
                      bool okay = aln->Align(seq.c_str(), ref_seq.c_str(),
                              ref_seq.length(), filter, &aln_result, 15);
                      if(!okay) continue;
                      uint16_t cur_score = aln_result.sw_score;
                      if((int)cur_score > best_score)
                      {
                          if(aln_result.mismatches == 0 && aln_result.cigar.size() == 1)
                          {
                            done = 1;
                          }
                          best_score = cur_score;
                          best_alignment = aln_result;
                          best_offset = ref_left_pos;
                          bestStrand = iter;
                      }
                  }    
              }
            }
            if(best_score > -1)
            {
                string refName = "*";
                size_t bestEnd = 0;
                size_t lastEnd = 0;
                map<size_t, string>::iterator it;

                for ( it = sapling->chrEnds.begin(); it != sapling->chrEnds.end(); it++ )
                {
                  size_t endCoord = it->first;
                  if(endCoord > best_alignment.ref_begin + best_offset && (bestEnd == 0 || endCoord < bestEnd))
                  {
                    bestEnd = endCoord;
                    refName = it->second;
                  }
                  if(endCoord <= best_alignment.ref_begin + best_offset && (lastEnd == 0 || endCoord > lastEnd))
                  {
                    lastEnd = endCoord;
                  }
                }
                best_alignment.ref_begin += best_offset - lastEnd;
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

int main(int argc, char** argv)
{
    if(argc < 4)
    {
        usage();
        return 1;
    }
    parseArgs(argc, argv);
    SaplingAligner* sa = new SaplingAligner(queryFn, refFn);
    sa->align_all_reads(outFn, argc, argv);
    delete sa;
    return 0;
}
