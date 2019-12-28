/*
 * SAPLING: Suffix Array Piecewise Linear INdex for Genomics
 */
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

#include "sa.h"
#include "util.h"

FILE* errorFile;

struct Sapling {

    string reference;
    int alpha = 2;
    int k = 21;
    int buckets = 18;
    double mostThreshold = 0.95;
    vector<size_t> sa; //sa[i] is the location in the suffix array where character i in reference appears
    vector<size_t> rev; // the inverse of sa sa[rev[i]] = i for all i
    size_t n;
    vector<int> overs, unders;
    int maxOver, maxUnder, meanError;
    int mostOver, mostUnder;
    size_t perfectPredictions = 0;

    map<string, int> chrEnds; // Map from chromosome to last position in genome which is part of that chromosome (zero-indexed)
    
    int vals[256];

    long long* xlist;
    long long* ylist;

    SuffixArray lsa;
    
    /*
     * Hashes the first k characters of a string into a 2k-bit integer
     */
    long long kmerize(const string& s)
    {
      long long kmer = 0;
      for(int i = 0; i<k; i++) kmer = (kmer << alpha) | vals[(size_t)s[i]];
      return kmer;
    }
    
    /*
     * Queries the piecewise linear index for a given kmer value.
     * Returns the predicted suffix array position
     */
    size_t queryPiecewiseLinear(long long x)
    {
      size_t bucket = (x >> (alpha*k - buckets));
      long long xlo = xlist[bucket];
      long long xhi = xlist[bucket+1];
      long long ylo = ylist[bucket];
      long long yhi = ylist[bucket+1];
      if(xlo == xhi) return ylo;
      long long predict = (long long)(.5 + ylo + (yhi - ylo) * (x - xlo) * 1. / (xhi - xlo));
      return (size_t)predict;
    }
    
    /*
     * Gets the lcp of a query string s with a given length 
     * and the reference starting at a particular position
     */
    size_t getLcp(size_t idx, string& s, int start, int length)
    {
      int i = start;
      for(; i<length && idx+i < n; i++) if(s[i] != reference[idx+i]) return i;
      return i; // whole thing matches
    }
  
    /*
     * Performs binary search in a section of the suffix array for suffixes starting with a given string
     *
     * s is the string being queried
     * lo is the start of the range
     * hi is the end of the range
     * loLcp is the LCP of s with the suffix at position lo
     * hiLcp is the LCP of s with the suffix at position hi
     * length is the number of characters in s
     * Returns a suffix array index such that the suffix there matches s if s occurs in the reference.
     */  
    long long binarySearch(string &s, size_t lo, size_t hi, size_t loLcp, size_t hiLcp, int length)
    {
      // Base case
      if(hi == lo + 2) return lo + 1;
  
      size_t mid = (lo + hi) >> 1;
      size_t idx = rev[mid];
      size_t nLcp = getLcp(idx, s, min(loLcp,  hiLcp), length);
      if(nLcp == s.length()) return (long long)mid;
      else if(lo + 1 >= hi) return (long long)-1;
      if(nLcp + idx == n || s[nLcp] > reference[idx+nLcp])
      {
        // suffix too small - search right half
        return binarySearch(s, mid, hi, nLcp, hiLcp, length);
      }
      else
      {
        // suffix too big - search left half
        return binarySearch(s, lo, mid, loLcp, nLcp, length);
      }
    }
    
    /*
     * Queries sapling for a given string s given its kmer value
     * The piecewise linear function is evaluated, and then a binary search is performed around the result
     */
    long long plQuery(string s, long kmer, size_t length)
    {
      //cout << s << endl;
      size_t predicted = queryPiecewiseLinear(kmer); // Predicted position in suffix array
      size_t idx = rev[predicted]; // Actual string position where we predict it to be
      size_t lcp = getLcp(idx, s, 0, length);
      if(lcp == length) return idx;
      size_t lo, hi;
      size_t loLcp = -1, hiLcp = -1;
      //cout << predicted << " " << idx << endl;
      if(lcp + idx == n || s[lcp] > reference[idx+lcp])
      {
        // Suffix is smaller then query - look farther right
        lo = predicted;
        hi = min(n-1, predicted+mostOver); // Over-prediction which the actual position is highly likely to not exceed
        size_t hiIdx = rev[hi]; // String index corresponding to over-prediction
        size_t oLcp = getLcp(hiIdx, s, 0, length); // LCP between over-prediction suffix and query
        if(oLcp == length) return hiIdx; // Over-prediction happened to be exactly right
        if(oLcp + hiIdx == n || s[oLcp] > reference[hiIdx+oLcp])
        {
          // bad case: over-prediction still not high enough
          lo = hi;
          loLcp = oLcp;
          hi = min(n-1, predicted + maxOver);
          hiIdx = rev[hi];
          oLcp = getLcp(hiIdx, s, 0, length);
          if(oLcp == length) return hiIdx;
          hiLcp = oLcp;
        }
        else
        {
          // correct position somewhere between original prediction and over-prediction
          loLcp = lcp;
          hiLcp = oLcp;
        }
      }
      else
      {
        // Suffix is bigger than query - look farther left
        lo = (size_t)max(0, (int)predicted-mostUnder);
        hi = predicted;
        size_t loIdx = rev[lo];
        size_t oLcp = getLcp(loIdx, s, 0, length); // LCP between under-prediction suffix and query
        if(oLcp == s.length()) return loIdx; // Under-prediction happened to be exactly right
        if(oLcp + loIdx == n || s[oLcp] > reference[loIdx+oLcp])
        {
          // correct position somewhere between original prediction and under-prediction
          hiLcp = lcp;
          loLcp = oLcp;
        }
        else
        {
          // bad case: under-prediction still not low enough
          hi = lo;
          hiLcp = oLcp;
          lo = (size_t)max(0, (int)predicted-maxUnder);
          loIdx = rev[lo];
          oLcp = getLcp(loIdx, s, 0, length);
          if(oLcp == s.length()) return loIdx;
          loLcp = oLcp;
        }
      }
      long long revPos = binarySearch(s, lo, hi, loLcp, hiLcp, length);
      if(revPos == -1) return -1;
      return rev[revPos];
    }

  size_t MAX_HITS = 32;

  /*
   * Counts the number of suffixes to the right of a position which have at least k characters in common with it
   * Used for getting all seeds for alignment
   */
  size_t countHitsRight(int sa_pos)
  {
    int lo = sa_pos, hi = sa_pos + MAX_HITS;
    if(hi > (int)(n - k)) hi = n - k + 1;
    while(hi > lo + 1)
    {
      int mid = (lo + hi)/2;
      if(lsa.queryLcpFromSAPos(sa_pos, mid))
        lo = mid;
      else
        hi = mid;
    }
    return lo - sa_pos;
   }

  /*
   * Counts the number of suffixes to the left of a position which have at least k characters in common with it
   * Used for getting all seeds for alignment
   */
  size_t countHitsLeft(int sa_pos)
  {
    int  lo = sa_pos - MAX_HITS, hi = sa_pos;
    if(lo < 0) lo = -1;
    while(hi > lo + 1)
    {
      int mid = (lo + hi)/2;
      if(lsa.queryLcpFromSAPos(mid, sa_pos))
        hi = mid;
      else
        lo = mid;
    }
    return sa_pos - hi;
  }
  
  /*
   * Calculates the prediction error given the prediction and true suffix array position
   * The true position can shift closer to the prediction when there are many occurrences of the same kmer
   */
  int getError(size_t y, size_t predict)
  {
    if(y < predict)
    {
      // Try increasing y to be closer to prediction
      long long lo = y, hi = (long long)predict + 1;
      while(lo < hi - 1)
      {
        size_t mid = (size_t)((lo + hi) / 2);
        if(lsa.queryLcpK(y, mid)) lo = (long long)mid;
        else hi = (long long)mid;
      }
      y = (size_t)lo;
      return (int)((long long)y - (long long)predict);
    }
    else if(y == predict) return 0;
    else 
    {
      long long lo = (long long)predict - 1, hi = y;
      while(lo < hi - 1)
      {
        size_t mid = (size_t)((lo + hi) / 2);
        if(lsa.queryLcpK(mid, y)) hi = (long long)mid;
        else lo = (long long)mid;
      }
      return (int)((long long)y - (long long)predict);
    }
  }

  /*
   * Calculates statistics about the prediction errors, such as the max and 95th percentile in each direction
   */
  void errorStats()
  {
    cout << "Computing error stats" << endl;
    maxUnder = 0;
    maxOver = 0;
    long tot = 0;
    size_t n = overs.size() + unders.size() + perfectPredictions;
    errorFile = fopen("errors.txt", "w");
    for(size_t i = 0; i<overs.size(); i++)
    {
      fprintf(errorFile, "%d\n", overs[i]);
      maxOver = max((int)overs[i], maxOver);
      tot += abs(overs[i]);
    }
    for(size_t i = 0; i<unders.size(); i++)
    {
      fprintf(errorFile, "%d\n", -unders[i]);
      maxUnder = max((int)unders[i], maxUnder);
      tot += abs(unders[i]);
    }
    if(maxUnder < 2) maxUnder = 2;  
    if(maxOver < 2) maxOver = 2;
    cout << "All overestimates within: " << maxOver << endl;
    cout << "All underestimates within: " << maxUnder << endl;
    cout << "Prefect predictions: " << perfectPredictions << endl;
    meanError = (int)(.5 + tot / n);
    cout << "Mean error: " << meanError << endl;
    sort(overs.begin(), overs.end());
    sort(unders.begin(), unders.end());
    mostOver = mostUnder = 0;
    if(overs.size() > 0) mostOver = overs[(size_t)(mostThreshold * overs.size())];
    if(unders.size() > 0) mostUnder = unders[(size_t)(mostThreshold * unders.size())];
    if(mostOver < 1) mostOver = 1;
    if(mostUnder < 1) mostUnder = 1;
    cout << mostThreshold << " of overestimates within: " << mostOver << endl;
    cout << mostThreshold << " of underestimates within: " << mostUnder << endl;
  }

  /*
   * Builds the piecewise linear index for a given string and suffix array
   */
  void buildPiecewiseLinear(string& s)
  {
    // Compute the maximum number of buckets < 5% of genome
    if(buckets == -1)
    {
      buckets = 1;
      while((size_t)(1L<<buckets) * 40 <= s.length()) buckets++;
    }
    printf("Buckets (log): %d\n", buckets);

    // Get hash of first kmer
    long long hash = kmerize(s.substr(0, k));

    // Initialize checkpoints
    xlist = new long long[(1L<<buckets)+1];
    ylist = new long long[(1L<<buckets)+1];
    for(size_t i = 0; i<(size_t)(1L<<buckets)+1; i++) xlist[i] = -1;
    for(size_t i = 0; i+k<=s.length(); i++)
    {
      // Grab hash as x value and then update the hash for next time
      long long x = hash;
      hash &= (1L << (2*(k-1))) - 1;
      hash <<= 2;
      if(i + k < s.length()) hash |= vals[(size_t)s[i+k]];

      // y is the suffix array position
      size_t y = sa[i];

      // See if this is a new checkpoint
      size_t bucket = (x >> (alpha*k - buckets));
      if(xlist[bucket] == -1 || xlist[bucket] > x)
      {
        xlist[bucket] = x;
        ylist[bucket] = y;
      }

      // Make sure the final checkpoint is the last point
      if(x > xlist[(1L<<buckets)])
      {
        xlist[(1L<<buckets)] = x;
        ylist[(1L<<buckets)] = y;
      }
    }

    // Fill in empty buckets
    if(xlist[0] == -1)
    {
      xlist[0] = 0;
      ylist[0] = 0;
    }
    for(size_t i = 1; i<(size_t)(1L<<buckets)+1; i++)
    {
      if(xlist[i] == -1)
      {
        xlist[i] = xlist[i-1];
        ylist[i] = ylist[i-1];
      }
    }

    // Make arrays for over and under predictions
    overs.resize(0);
    unders.resize(0);
    perfectPredictions = 0;

    // Initialize the hash to sweep through kmers again
    hash = kmerize(s.substr(0, k));
    for(size_t i = 0; i+k<=s.length(); i++)
    {
      // Get predicted suffix array position
      size_t predict = queryPiecewiseLinear(hash);

      // Update hash for next time
      hash &= (1L << (2*(k-1))) - 1;
      hash <<= 2;
      if(i + k < s.length()) hash |= vals[(size_t)s[i+k]];

      // Get actual suffix array position and figure out the error
      size_t y = sa[i];
      int val = getError(y, predict);

      // Update corresponding list of errors
      if(val > 0) overs.push_back(val);
      else if(val < 0)
      {
        unders.push_back(-val);
      }
      else perfectPredictions++;
    }

    // Compute maximum/average/etc. errors
    errorStats();
  }
  
  /*
   * Takes a FASTA filepath and builds Sapling from the contained genome
   */
  Sapling(string refFnString, string saFnString, string saplingFnString, int numBuckets)
  {
    for(int i = 0; i<256; i++) vals[i] = 0;
    vals['A'] = 0;
    vals['C'] = 1;
    vals['G'] = 2;
    vals['T'] = 3;

    buckets = numBuckets;
      
    ifstream input(refFnString);
    string cur;
    std::ostringstream out("");
    cout << "Reading reference genome" << endl;
    size_t charCount = 0;
    string curName = "";
    while (getline(input, cur))
    {
      if(cur[0] != '>')
      {
        for(size_t i = 0; i<cur.length(); i++)
        {
          if(cur[i] >= 'a' && cur[i] <= 'z') cur[i] += 'A' - 'a';
          if(!bad(cur[i]))
          {
            charCount++;
            out << cur[i];
          }
        }
      }
      else
      {
        if(curName.length() > 0)
        {
          chrEnds[curName] = charCount;
        }
        auto first_token = cur.substr(0, cur.find(' '));
        curName = first_token.substr(1);
      }
    }
    if(curName.length() > 0)
    {
      chrEnds[curName] = charCount;
    }
    reference = out.str();
      
    n = reference.length();

    const char *fn = saFnString.c_str();
    ifstream f(fn);
    if(f.good())
    {
      cout << "Reading suffix array from file" << endl;
      FILE *infile = fopen (fn, "rb");
      size_t size;
      size_t err = fread(&size, sizeof(size_t), 1, infile);
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
      lsa = SuffixArray();
          
      cout << "Constructing RMQ" << endl;
      lsa.lcp = lcp;
      lsa.krmq_init(k);
      lsa.inv = sa;
          
      cout << "Loaded suffix array of size " << sa.size() << endl;
    }
    else
    {
      cout << "Building suffix array" << endl;
      lsa = sa_init3(reference, alpha);
      cout << "Writing suffix array to file" << endl;
      sa = lsa.inv;
      FILE *outfile = fopen (fn, "wb");
      size_t size = sa.size();
      fwrite(&size, sizeof(size_t), 1, outfile);
      fwrite(&sa[0], sizeof(size_t), size, outfile);
      size = lsa.lcp.size();
      fwrite(&size, sizeof(size_t), 1, outfile);
      fwrite(&lsa.lcp[0], sizeof(size_t), size, outfile);
      cout << "Making LCP RMQ" << endl;
      lsa.krmq_init(k);
      cout << "Built suffix array of size " << sa.size() << endl;
    }

    vector<size_t>().swap(lsa.lcp);
    
    cout << "Initializing rev and sa" << endl;
    rev = vector<size_t>(n, 0);
    cout << "Filling rev and sa" << endl;
    for(size_t i = 0; i<n; i++) rev[sa[i]] = i;
    
    const char *saplingfn = saplingFnString.c_str();
    ifstream saplingf(saplingfn);
    if(saplingf.good())
    {
      cout << "Reading Sapling from file" << endl;
      FILE *infile = fopen (saplingfn, "rb");
      size_t xlistsize;
      int err = fread(&buckets, sizeof(int), 1, infile);
      if(buckets <= 30)
      {
        int xlsize;
        int err = fread(&xlsize, sizeof(int), 1, infile);
        if(err != 1)
        {
          cerr << "Error reading sapling data structure fom file" << endl;
        }
        xlistsize = (size_t)xlsize;
      }
      else
      {
        int err = fread(&xlistsize, sizeof(size_t), 1, infile);
        if(err != 1)
        {
          cerr << "Error reading sapling data structure from file" << endl;
        }
      }
      xlist = new long long[xlistsize];
      ylist = new long long[xlistsize];
      err = fread(&xlist[0], sizeof(long long), xlistsize, infile);
      err = fread(&ylist[0], sizeof(long long), xlistsize, infile);
      err = fread(&maxOver, sizeof(int), 1, infile);
      err = fread(&maxUnder, sizeof(int), 1, infile);
      err = fread(&meanError, sizeof(int), 1, infile);
      err = fread(&mostOver, sizeof(int), 1, infile);
      err = fread(&mostUnder, sizeof(int), 1, infile);
      if(err != 1)
      {
        cerr << "Error reading sapling data structure from file" << endl;
      }
    }
    else
    {
      cout << "Building Sapling" << endl;
      buildPiecewiseLinear(reference);
      cout << "Writing Sapling to file" << endl;
      FILE *outfile = fopen (saplingfn, "wb");
      size_t xlistsize = (1L<<buckets)+1;
      fwrite(&buckets, sizeof(int), 1, outfile);
      if(buckets <= 30)
      {
        fwrite(&xlistsize, sizeof(int), 1, outfile);
      }
      else
      {
        int xlsize = (1<<buckets)+1;
        fwrite(&xlsize, sizeof(int), 1, outfile);
      }
      fwrite(&xlist[0], sizeof(long long), xlistsize, outfile);
      fwrite(&ylist[0], sizeof(long long), xlistsize, outfile);
      fwrite(&maxOver, sizeof(int), 1, outfile);
      fwrite(&maxUnder, sizeof(int), 1, outfile);
      fwrite(&meanError, sizeof(int), 1, outfile);
      fwrite(&mostOver, sizeof(int), 1, outfile);
      fwrite(&mostUnder, sizeof(int), 1, outfile);
    }
  }
};
