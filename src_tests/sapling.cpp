/*
 * SAPLING: Suffix Array Piecewise Linear INdex for Genomics
 * Note: compile with -O2 -std=c++11 for timing function to work
 */
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>
#include "sa.h"
#include "util.h"
 
using namespace std::chrono;

string reference;
int replace = 1; // whether or not to replace all non-base characters with A's
int alpha = 2;
int k = 21;
int buckets = 15; // TODO handle cases when some buckets are empty
double mostThreshold = 0.95;
vector<size_t> sa; //sa[i] is the location in the suffix array where character i in reference appears
vector<size_t> rev;
size_t n;
vector<long long> errors;
vector<int> overs, unders;
int vals[256];
int maxOver, maxUnder, meanError;
int mostOver, mostUnder;

long long* xlist;
long long* ylist;

SuffixArray lsa;

string decode(long long val)
{
  string res = "";
  for(int i = 0; i<k; i++)
  {
    long long cur = val%4;
    string c = "T";
    if(cur == 0) c = "A";
    else if(cur == 1) c = "C";
    else if(cur == 2) c = "G";  
    res = c + res;
    val /= 4;
  }
  return res;
}

long long kmerize(const string& s)
{
    long long kmer = 0;
	for(int i = 0; i<k; i++) kmer = (kmer << alpha) | vals[s[i]];
	return kmer;
}

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

size_t getLcp(size_t idx, string& s, int start, int length)
{
	int i = start;
	for(; i<length && idx+i < n; i++) if(s[i] != reference[idx+i]) return i;
	return i; // whole thing matches
}

size_t binarySearch(string &s, size_t lo, size_t hi, size_t loLcp, size_t hiLcp, int length)
{
	// Base case
	if(hi == lo + 2) return lo + 1;
	
	size_t mid = (lo + hi) >> 1;
	size_t idx = rev[mid];
	size_t nLcp = getLcp(idx, s, min(loLcp,  hiLcp), length);
	if(nLcp == s.length()) return mid;
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

size_t plQuery(string &s, long kmer, int length)
{
	size_t predicted = queryPiecewiseLinear(kmer); // Predicted position in suffix array
	size_t idx = rev[predicted]; // Actual string position where we predict it to be
	size_t lcp = getLcp(idx, s, 0, length);
	if(lcp == length) return idx;
	size_t lo, hi;
	size_t loLcp = -1, hiLcp = -1;
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
	return rev[binarySearch(s, lo, hi, loLcp, hiLcp, length)];
}

int getError(size_t y, size_t predict)
{
	if(y < predict)
	{
		size_t lo = y, hi = predict+1;
		while(lo < hi - 1)
		{
			size_t mid = (lo+hi)/2;
			if(lsa.queryLcpK(mid, y))
			{
				lo = mid;
			}
			else hi = mid;
		}
		return (int)lo - (int)predict;
	}
	else if(y == predict) return 0;
	else 
	{
		size_t lo = predict - 1, hi = y;
		while(lo < hi - 1)
		{
			size_t mid = (lo+hi)/2;
			if(lsa.queryLcpK(mid, y))
			{
				hi = mid;
			}
			else lo = mid;
		}
		return (int)hi - (int)predict;
	}
}

void errorStats()
{
  cout << "Computing error stats" << endl;
	maxUnder = 0;
	maxOver = 0;
	long tot = 0;
	size_t n = overs.size() + unders.size();
	for(size_t i = 0; i<overs.size(); i++)
	{
    //cout << "over: " << overs[i] << endl;
		maxOver = max((int)overs[i], maxOver);
		tot += abs(overs[i]);
	}
  for(size_t i = 0; i<unders.size(); i++)
	{
    //cout << "under: " << unders[i] << endl;
		maxUnder = max((int)unders[i], maxUnder);
		tot += abs(unders[i]);
	}
	cout << "All overestimates within: " << maxOver << endl;
	cout << "All underestimates within: " << maxUnder << endl;
	meanError = (int)(.5 + tot / n);
	cout << "Mean error: " << meanError << endl;
	sort(overs.begin(), overs.end());
	sort(unders.begin(), unders.end());
	mostOver = overs[(size_t)(mostThreshold * overs.size())];
	mostUnder = unders[(size_t)(mostThreshold * unders.size())];
	cout << mostThreshold << " of overestimates within: " << mostOver << endl;
	cout << mostThreshold << " of underestimates within: " << mostUnder << endl;
}

void buildPiecewiseLinear(string& s, vector<size_t> sa)
{
	vector<long long> xs;
	vector<size_t> ys;
	long long hash = kmerize(s.substr(0, k));
	for(size_t i = 0; i+k<=s.length(); i++)
	{
		xs.push_back(hash);
		ys.push_back(sa[i]);
		hash &= (1L << (2*(k-1))) - 1;
		hash <<= 2;
		if(i + k < s.length()) hash |= vals[s[i+k]];
	}
	xlist = new long long[(1L<<buckets)+1];
	ylist = new long long[(1L<<buckets)+1];
	for(int i = 0; i<(1L<<buckets)+1; i++) xlist[i] = -1;
	for(size_t i = 0; i<xs.size(); i++)
	{
		long long x = xs[i];
		size_t y = ys[i];
		size_t bucket = 0;
		bucket = (x >> (alpha*k - buckets));
		if(xlist[bucket] == -1 || xlist[bucket] > x)
		{
			xlist[bucket] = x;
			ylist[bucket] = y;
		}
		if(x > xlist[(1L<<buckets)])
		{
			xlist[(1L<<buckets)] = x;
			ylist[(1L<<buckets)] = y;
		}
	}
	if(xlist[0] == -1)
	{
	    xlist[0] = 0;
	    ylist[0] = 0;
	}
	for(int i = 1; i<(1L<<buckets)+1; i++)
	{
	    if(xlist[i] == -1)
	    {
	        xlist[i] = xlist[i-1];
	        ylist[i] = ylist[i-1];
	    }
	}
	//errors.resize(xs.size());
	overs.resize(0);
	unders.resize(0);
  printf("(%ld, %ld) to (%ld, %ld)\n", xlist[0], ylist[0], xlist[1], ylist[1]);
  printf("(%ld, %ld) to (%ld, %ld)\n", xlist[(1L<<buckets)+1-2], ylist[(1L<<buckets)+1-2], xlist[(1L<<buckets)+1-1], ylist[(1L<<buckets)+1-1]);
  FILE *tmpFout1 = fopen("errorsfirst.txt", "w"), *tmpFout2 = fopen("errorslast.txt", "w"), *tmpFout3 = fopen("allerrors.txt", "w");
	for(size_t i = 0; i<xs.size(); i++)
	{
		size_t predict = queryPiecewiseLinear(xs[i]);
		size_t y = ys[i];

		int curError = getError(y, predict);
		if(curError > 0) overs.push_back(curError);
		else unders.push_back(-curError);

    fprintf(tmpFout3, "%ld %ld %ld %ld\n", xs[i], ys[i], predict, curError);
    
    if(xs[i] > xlist[1] && xs[i] < xlist[(1L<<buckets)+1-2])
    {
      //fprintf(tmpFout3, "%ld %ld %ld %ld\n", xs[i], ys[i], predict, errors[i]);
    }
    else if(xs[i] <= xlist[1])
    {
      string d = decode(xs[i]);
      //cout << i << " " << xs[i] << " " << ys[i] << " " << reference.substr(i, k) << " " << d << endl;

      //fprintf(tmpFout1, "%ld %ld %ld\n", xs[i], ys[i], predict);
    }
    else
    {
      //fprintf(tmpFout2, "%ld %ld %ld\n", xs[i], ys[i], predict);
    }
	}
  fclose(tmpFout1);
  fclose(tmpFout2);
  fclose(tmpFout3);
	errorStats();
}

int main(int argc, char **argv)
{
    if(argc != 5)
    {
        cout << "Usage: " << argv[0] << " <genome> " << " <suffix array file> " << " <sapling file> " << " <log number buckets> " << endl;
        return 0;
    }
    buckets = atoi(argv[4]);
    for(int i = 0; i<256; i++) vals[i] = 0;
    vals['A'] = (1<<alpha)-4;
    vals['C'] = (1<<alpha)-3;
    vals['G'] = (1<<alpha)-2;
    vals['T'] = (1<<alpha)-1;
    vals['N'] = (1<<alpha)-5;
    ifstream input(argv[1]);
    string cur;
    std::ostringstream out("");
    cout << "Reading reference genome" << endl;
    while (getline(input, cur))
    {
        if(cur[0] != '>')
        {
            for(int i = 0; i<cur.length(); i++)
            {
                if(cur[i] >= 'a' && cur[i] <= 'z') cur[i] += 'A' - 'a';
                if(!bad(cur[i]))
                   out << cur[i];
            }
        }
    }
    reference = out.str();
    n = reference.length();
    cout << n << endl;
    
    vector<size_t> sa;
    
    string fnString = argv[2];
    const char *fn = fnString.c_str();
    ifstream f(fn);
    if(f.good())
    {
        cout << "Reading suffix array from file" << endl;
        FILE *infile = fopen (fn, "rb");
        size_t size;
        fread(&size, sizeof(size_t), 1, infile);
        sa.resize(size);
        fread(&sa[0], sizeof(size_t), size, infile);
        
        vector<size_t> lcp;
        fread(&size, sizeof(size_t), 1, infile);
        lcp.resize(size);
        fread(&lcp[0], sizeof(size_t), size, infile);
        lsa = SuffixArray();
        
        cout << "Constructing RMQ" << endl;
        lsa.krmq = KRMQ(lcp, k);
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
        lsa.krmq = KRMQ(lsa.lcp, k);
        cout << "Built suffix array of size " << sa.size() << endl;
    }
    vector<size_t>().swap(lsa.lcp);
    
    cout << "Initializing rev and sa" << endl;
    rev = vector<size_t>(n, 0);
    cout << "Filling rev and sa" << endl;
    for(size_t i = 0; i<n; i++) rev[sa[i]] = i;
    
    string saplingfnString = argv[3];
    const char *saplingfn = saplingfnString.c_str();
    ifstream saplingf(saplingfn);
    if(saplingf.good())
    {
        cout << "Reading Sapling from file" << endl;
        FILE *infile = fopen (saplingfn, "rb");
        size_t xlistsize;
        if(buckets <= 30)
        {
            int xlsize;
            fread(&xlsize, sizeof(int), 1, infile);
            xlistsize = (size_t)xlsize;
        }
        else
        {
            fread(&xlistsize, sizeof(size_t), 1, infile);
        }
        xlist = new long long[xlistsize];
        ylist = new long long[xlistsize];
        fread(&xlist[0], sizeof(long long), xlistsize, infile);
        fread(&ylist[0], sizeof(long long), xlistsize, infile);
        fread(&maxOver, sizeof(int), 1, infile);
        fread(&maxUnder, sizeof(int), 1, infile);
        fread(&meanError, sizeof(int), 1, infile);
        fread(&mostOver, sizeof(int), 1, infile);
        fread(&mostUnder, sizeof(int), 1, infile);
    }
    else
    {
        cout << "Building Sapling" << endl;
        buildPiecewiseLinear(reference, sa);
        cout << "Writing Sapling to file" << endl;
        FILE *outfile = fopen (saplingfn, "wb");
        size_t xlistsize = (1L<<buckets)+1;
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
    
    cout << "Testing Sapling" << endl;
    int numQueries = 5000000;
    vector<string> queries(numQueries, "");
    vector<long long> kmers = vector<long long>(numQueries, 0);
    for(int i = 0; i<numQueries; i++)
    {
    	size_t idx = rand() % (n - k);
    	queries[i] = reference.substr(idx, k);
    	kmers[i] = kmerize(queries[i]);
    }
    
    FILE *outfile = fopen ("queries.out", "w");
    for(int i = 0; i < numQueries; i++)
    {
        fprintf(outfile, "@read%d\n", i+1);
        fprintf(outfile, "%s\n", queries[i].c_str());
        fprintf(outfile, "+\n");
        
        for(int j = 0; j<k; j++) fprintf(outfile, "9");
        fprintf(outfile, "\n");
    }
    
    cout << "Constructed queries" << endl;
    // Run piece-wise linear test
    vector<size_t> plAnswers(numQueries, 0);
    auto start = std::chrono::system_clock::now();
    for(int i = 0; i<numQueries; i++)
    {
        plAnswers[i] = plQuery(queries[i], kmers[i], queries[i].length());
    }
    
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "Piecewise linear time: " << elapsed_seconds.count() << endl;
    
    // Check the answers
    int countCorrect = 0;
    for(int i = 0; i<numQueries; i++)
    {
        if(plAnswers[i] + k <= n && queries[i] == reference.substr(plAnswers[i], k))
        {
            countCorrect++;
        }
    }
    cout <<"Piecewise linear correctness: " << countCorrect << " out of " << numQueries << endl;
}
