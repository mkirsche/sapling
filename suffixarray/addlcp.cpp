using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <vector>

size_t size;
vector<size_t> inv, idx;
vector<char> str;

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

vector<size_t> getLCP() 
{
  inv.resize(size);
  for(size_t i = 0; i<size; i++)
  {
    inv[idx[i]] = i;
  }
  vector<size_t> lcp(size - 1, 0);
  size_t curr = 0;
  for (size_t i = 0; i < size; i++) 
  {
   // printf("%d %d %c\n", i, inv[i], str[i]);
    size_t k = inv[i];
    if (k < size - 1)
    {
      size_t j = idx[k + 1];
      //printf("j %d %c\n", j, str[j]);
      while (i + curr < size && j + curr < size && str[i + curr] == str[j + curr])
      {
        curr++;
      }
      //printf("curr %d\n", curr);
      lcp[k] = curr;
    }
    if (curr > 0)
    {
      curr--;
    }
  }
  printf("lcp done\n");
  return lcp;
}

int main(int argc, const char *argv[])
{
  printf("reference: %s\n", argv[1]);
  size = filesize(argv[1]);
  str.resize(size);
  printf("n: %ld\n", size);
  FILE *fp = fopen(argv[1], "rb");
  fread(&str[0], 1, size, fp);

  idx.resize(size);

  fp = fopen(argv[2], "rb");
  fread(&idx[0], sizeof(size_t), size, fp);

  fclose(fp);
  vector<size_t> lcp = getLCP();

  FILE *outfile = fopen (argv[3], "wb");
  fwrite(&size, sizeof(size_t), 1, outfile);
  fwrite(&inv[0], sizeof(size_t), size, outfile);
  size_t lcpSize = lcp.size();
  fwrite(&lcpSize, sizeof(size_t), 1, outfile);
  fwrite(&lcp[0], sizeof(size_t), lcpSize, outfile);

  return 0;
}
