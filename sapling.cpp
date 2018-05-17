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
vector<int> sa; //sa[i] is the location in the suffix array where character i in reference appears
vector<int> rev; // the inverse of sa sa[rev[i]] = i for all i
int n;
vector<int> errors;
vector<int> overs, unders;
int vals[256];
int maxOver, maxUnder, meanError;
int mostOver, mostUnder;

long long* xlist;
long long* ylist;

SuffixArray lsa;

long long kmerize(const string& s)
{
    long long kmer = 0;
	for(int i = 0; i<k; i++) kmer = (kmer << alpha) | vals[s[i]];
	return kmer;
}

int queryPiecewiseLinear(long long x)
{
	int bucket = (int)(x >> (alpha*k - buckets));
	long long xlo = xlist[bucket];
	long long xhi = xlist[bucket+1];
	long long ylo = ylist[bucket];
	long long yhi = ylist[bucket+1];
	long long predict = (long long)(.5 + ylo + (yhi - ylo) * (x - xlo) * 1. / (xhi - xlo));
	return (int)predict;
}

int getLcp(int idx, string& s, int start, int length)
{
	int i = start;
	for(; i<length && idx+i < n; i++) if(s[i] != reference[idx+i]) return i;
	return i; // whole thing matches
}

int binarySearch(string &s, int lo, int hi, int loLcp, int hiLcp, int length)
{
	// Base case
	if(hi == lo + 2) return lo + 1;
	
	int mid = (lo + hi) >> 1;
	int idx = rev[mid];
	int nLcp = getLcp(idx, s, min(loLcp,  hiLcp), length);
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

int plQuery(string &s, long kmer, int length)
{
	int predicted = queryPiecewiseLinear(kmer); // Predicted position in suffix array
	int idx = rev[predicted]; // Actual string position where we predict it to be
	int lcp = getLcp(idx, s, 0, length);
	if(lcp == length) return idx;
	int lo, hi;
	int loLcp = -1, hiLcp = -1;
	if(lcp + idx == n || s[lcp] > reference[idx+lcp])
	{
		// Suffix is smaller then query - look farther right
		lo = predicted;
		hi = min(n-1, predicted+mostOver); // Over-prediction which the actual position is highly likely to not exceed
		int hiIdx = rev[hi]; // String index corresponding to over-prediction
		int oLcp = getLcp(hiIdx, s, 0, length); // LCP between over-prediction suffix and query
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
		lo = max(0, predicted-mostUnder);
		hi = predicted;
		int loIdx = rev[lo];
		int oLcp = getLcp(loIdx, s, 0, length); // LCP between under-prediction suffix and query
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
			lo = max(0, predicted - maxUnder);
			loIdx = rev[lo];
			oLcp = getLcp(loIdx, s, 0, length);
			if(oLcp == s.length()) return loIdx;
			loLcp = oLcp;
		}
	}
	return rev[binarySearch(s, lo, hi, loLcp, hiLcp, length)];
}

int getError(int y, int predict)
{
	if(y < predict)
	{
		int lo = y, hi = predict+1;
		while(lo < hi - 1)
		{
			int mid = (lo+hi)/2;
			if(lsa.queryLcp(mid, y) >= k)
			{
				lo = mid;
			}
			else hi = mid;
		}
		return lo - predict;
	}
	else if(y == predict) return 0;
	else 
	{
		int lo = predict - 1, hi = y;
		while(lo < hi - 1)
		{
			int mid = (lo+hi)/2;
			if(lsa.queryLcp(mid, y) >= k)
			{
				hi = mid;
			}
			else lo = mid;
		}
		return hi - predict;
	}
}

void errorStats()
{
    cout << "Computing error stats" << endl;
	maxUnder = 0;
	maxOver = 0;
	long tot = 0;
	int n = overs.size() + unders.size();
	for(int i = 0; i<n; i++)
	{
		maxUnder = max(-errors[i], maxUnder);
		maxOver = max(errors[i], maxOver);
		tot += abs(errors[i]);
	}
	cout << "All overestimates within: " << maxOver << endl;
	cout << "All underestimates within: " << maxUnder << endl;
	meanError = (int)(.5 + tot / n);
	cout << "Mean error: " << meanError << endl;
	sort(overs.begin(), overs.end());
	sort(unders.begin(), unders.end());
	mostOver = overs[(int)(mostThreshold * overs.size())];
	mostUnder = unders[(int)(mostThreshold * unders.size())];
	cout << mostThreshold << " of overestimates within: " << mostOver << endl;
	cout << mostThreshold << " of underestimates within: " << mostUnder << endl;
}

void buildPiecewiseLinear(string& s, vector<int> sa)
{
	vector<long long> xs;
	vector<int> ys;
	long long hash = kmerize(s.substr(0, k));
	for(int i = 0; i+k<=s.length(); i++)
	{
		xs.push_back(hash);
		ys.push_back(sa[i]);
		hash &= (1L << (2*(k-1))) - 1;
		hash <<= 2;
		if(i + k < s.length()) hash |= vals[s[i+k]];
	}
	xlist = new long long[(1<<buckets)+1];
	ylist = new long long[(1<<buckets)+1];
	for(int i = 0; i<(1<<buckets)+1; i++) xlist[i] = -1;
	for(int i = 0; i<xs.size(); i++)
	{
		long long x = xs[i];
		int y = ys[i];
		int bucket = 0;
		bucket = (int)(x >> (alpha*k - buckets));
		if(xlist[bucket] == -1 || xlist[bucket] > x)
		{
			xlist[bucket] = x;
			ylist[bucket] = y;
		}
		if(x > xlist[(1<<buckets)])
		{
			xlist[(1<<buckets)] = x;
			ylist[(1<<buckets)] = y;
		}
	}
	errors.resize(xs.size());
	overs.resize(0);
	unders.resize(0);
	for(int i = 0; i<xs.size(); i++)
	{
		int predict = queryPiecewiseLinear(xs[i]);
		int y = ys[i];
		errors[i] = getError(y, predict);
		if(errors[i] > 0) overs.push_back(errors[i]);
		else unders.push_back(-errors[i]);
	}
	errorStats();
}

int main()
{
	for(int i = 0; i<256; i++) vals[i] = 0;
	vals['A'] = (1<<alpha)-4;
	vals['C'] = (1<<alpha)-3;
	vals['G'] = (1<<alpha)-2;
	vals['T'] = (1<<alpha)-1;
	vals['N'] = (1<<alpha)-5;
	ifstream input("chr22.fa");
	string s;
	string cur;
	std::ostringstream out("");
	cout << "Reading reference genome" << endl;
	while (getline(input, cur))
	{
		if(cur[0] != '>') out << cur;
	}
	s = out.str();
	cout << "Removing non-base characters" << endl;
    s = rep(s);
    reference = s;
    n = reference.length();
    cout << s.length() << endl;
	
    vector<int> x;
	
    const char *fn = "sa_chr22.txt";
    ifstream f(fn);
    if(f.good())
    {
        cout << "Reading suffix array from file" << endl;
        FILE *infile = fopen (fn, "rb");
        size_t size;
        fread( &size, sizeof( size_t ), 1, infile);
        x.resize(size);
        fread(&x[0], sizeof(int), size, infile);
        
        vector<int> lcp;
        fread( &size, sizeof( size_t ), 1, infile);
        lcp.resize(size);
        fread(&lcp[0], sizeof(int), size, infile);
        lsa = SuffixArray();
        lsa.rmq = RMQ(lcp);
        lsa.inv = x;
    }
    else
    {
        cout << "Building suffix array" << endl;
	    lsa = sa_init3(s, alpha);
	    cout << "Writing suffix array to file" << endl;
	    x = lsa.inv;
	    FILE *outfile = fopen (fn, "wb");
	    size_t size = x.size();
	    fwrite( &size, sizeof( size_t ), 1, outfile);
	    fwrite( &x[0], sizeof(int), size, outfile);
	    
	    size = lsa.lcp.size();
	    fwrite( &size, sizeof( size_t ), 1, outfile);
	    fwrite( &lsa.lcp[0], sizeof(int), size, outfile);
	}
	cout << "Built suffix array of size " << x.size() << endl;
	sa = vector<int>(n, 0); 
	rev = vector<int>(n, 0);
	cout << "Initialized rev and sa" << endl;
	for(int i = 0; i<n; i++) rev[sa[i] = x[i]] = i;
	cout << "Building Sapling" << endl;
	buildPiecewiseLinear(reference, x);
	cout << "Testing Sapling" << endl;
	int numQueries = 5000000;
	vector<string> queries(numQueries, "");
	vector<long long> kmers = vector<long long>(numQueries, 0);
	for(int i = 0; i<numQueries; i++)
	{
		int idx = rand() % (n - k);
		queries[i] = s.substr(idx, k);
		kmers[i] = kmerize(queries[i]);
	}
	cout << "Constructed queries" << endl;
	// Run piece-wise linear test
	vector<int> plAnswers(numQueries, 0);
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
		if(plAnswers[i] + k <= n && queries[i] == s.substr(plAnswers[i], k))
		{
			countCorrect++;
		}
	}
	cout <<"Piecewise linear correctness: " << countCorrect << " out of " << numQueries << endl;
}
