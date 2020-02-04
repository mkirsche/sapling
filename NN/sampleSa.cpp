using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int vals[256];

vector<char> genome;

long long kmerize(size_t start, int k)
{
  long long kmer = 0;
	for(size_t i = start; i<(size_t)(k) + start; i++)
  {
    kmer = (kmer << 2) | vals[genome[i]];
  }
	return kmer;
}

string genome_substr(size_t start, int k)
{
  string kmer = "";
	for(size_t i = start; i<(size_t)(k) + start; i++)
  {
    kmer += genome[i];
  }
	return kmer;
}

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

int k = 21;

int main(int argc, const char *argv[])
{
  for(int i = 0; i<256; i++) vals[i] = 0;
  vals['A'] = 0;
  vals['C'] = 1;
  vals['G'] = 2;
  vals['T'] = 3;

  size_t size = filesize(argv[1]) / 8;
  vector<size_t> idx, inv;
  idx.resize(size);
  inv.resize(size);

  FILE *fp = fopen(argv[1], "rb");
  fread(&idx[0], sizeof(size_t), size, fp);
  fclose(fp);
  for(size_t i = 0; i < size; i++) inv[idx[i]] = i;

  genome.resize(size);
  fp = fopen(argv[2], "rb");
  fread(&genome[0], sizeof(char), size, fp);

  size_t sample = 1;

  fclose(fp);

  for(size_t i = 0; i<size - k + 1; i+=sample)
  {
    cout << inv[i] << " " << kmerize(i, k) << endl;
  }

  return 0;
}
