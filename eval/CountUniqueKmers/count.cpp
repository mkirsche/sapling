/*
 * Uses the LCP array to count the number of different k-mers and how many of them are unique
 */
using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, const char *argv[])
{
  if(argc != 2)
  {
    cout << "Usage: " << argv[0] << " <suffix array file>" << endl;
    return 1;
  }
  string saFnString = argv[1];
  cerr << "Suffix Array File: " << saFnString << endl;
  const char *fn = saFnString.c_str();
  ifstream f(fn);
  if(f.good())
  {
    cerr << "Reading suffix array from file" << endl;
    FILE *infile = fopen (fn, "rb");
    size_t size;
    size_t err = fread(&size, sizeof(size_t), 1, infile);
    vector<size_t> sa;
    sa.resize(size);
    err = fread(&sa[0], sizeof(size_t), size, infile);
    
    vector<size_t> lcp;
    err = fread(&size, sizeof(size_t), 1, infile);
    lcp.resize(size);
    err = fread(&lcp[0], sizeof(size_t), size, infile);
    if(err == 0)
    {
      cerr << "Error reading suffix array from file" << endl;
    }
    cerr << "Loaded suffix array of size " << sa.size() << endl;
    vector<size_t> counts(1001, 0);
    vector<size_t> counts2(counts.size(), 0);
    for(size_t i = 0; i<size-1; i++)
    {
      //cout << i << " " << lcp[i] << endl;
      if(lcp[i] >= counts.size())
      {
        //cout << lcp[i] << endl;
        lcp[i] = counts.size()-1;
      }
      if(i < size - 2 && lcp[i+1] >= counts.size())
      {
        //cout << lcp[i] << endl;
        lcp[i+1] = counts.size()-1;
      }
      if(i < size - 2)
      {
        counts2[max(lcp[i], lcp[i+1])]++;
      }
      counts[lcp[i]]++;
    }
    cout << "counting done" << endl;
    size_t cumulativeSum = 0;
    size_t cumulativeSum2 = 0;
    for(size_t i = 1; i<counts.size(); i++)
    {
      if(counts2[i-1] == 0) counts2[i-1]++;
      if(counts[i-1] == 0) counts[i-1]++;
      cumulativeSum += counts[i-1];
      cumulativeSum2 += counts2[i-1];
      cout << i << " " << (cumulativeSum+1) << " " << (cumulativeSum2-1) << " " << (size-i+1) << endl;
      cumulativeSum--;
      cumulativeSum2--;
    }
  }
  else
  {
    cout << "Suffix Array File Not Found" << endl;
  }
  return 0;
}
