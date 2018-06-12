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
vector<size_t> rev; // the inverse of sa sa[rev[i]] = i for all i
size_t n;
vector<int> errors;
vector<size_t> overs, unders;
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

size_t queryPiecewiseLinear(long long x)
{
	int bucket = (int)(x >> (alpha*k - buckets));
	long long xlo = xlist[bucket];
	long long xhi = xlist[bucket+1];
	long long ylo = ylist[bucket];
	long long yhi = ylist[bucket+1];
	long long predict = (long long)(.5 + ylo + (yhi - ylo) * (x - xlo) * 1. / (xhi - xlo));
	return (size_t)predict;
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
	for(size_t i = 0; i<n; i++)
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
	mostOver = overs[(size_t)(mostThreshold * overs.size())];
	mostUnder = unders[(size_t)(mostThreshold * unders.size())];
	cout << mostThreshold << " of overestimates within: " << mostOver << endl;
	cout << mostThreshold << " of underestimates within: " << mostUnder << endl;
}

void buildPiecewiseLinear(string& s, vector<size_t> sa, FILE *outfile, int printbuckets)
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
	xlist = new long long[(1<<buckets)+1];
	ylist = new long long[(1<<buckets)+1];
	for(int i = 0; i<(1<<buckets)+1; i++) xlist[i] = -1;
	for(size_t i = 0; i<xs.size(); i++)
	{
		long long x = xs[i];
		size_t y = ys[i];
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
	fprintf(outfile, "%s\t%s\t%s\t%s\n", "Error", "Kmer", "saPos", "Bucket");
	for(size_t i = 0; i<xs.size(); i++)
	{
		size_t predict = queryPiecewiseLinear(xs[i]);
		size_t y = ys[i];
		errors[i] = getError(y, predict);
		int bucket = (int)(xs[i] >> (alpha*k - buckets));
		if(printbuckets == -1 || bucket < printbuckets)
		{
		    fprintf(outfile, "%d\t%lld\t%zu\t%d\n", errors[i], xs[i], ys[i], bucket);
		}
		if(errors[i] > 0) overs.push_back(errors[i]);
		else unders.push_back(-errors[i]);
	}
	errorStats();
}

int main(int argc, char **argv)
{
    if(argc != 6)
    {
        cout << "Usage: " << argv[0] << " <genome> " << " <suffix array file> " << " <errors file> " 
            << " <log number buckets> " << " <buckets to print> " << endl;
        return 0;
    }
    buckets = atoi(argv[4]);
    int printbuckets = atoi(argv[5]);
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
    
    cout << "Initializing rev and sa" << endl;
    rev = vector<size_t>(n, 0);
    cout << "Filling rev and sa" << endl;
    for(size_t i = 0; i<n; i++) rev[sa[i]] = i;
    
    cout << "Building Sapling and writing errors to file" << endl;
    string errorsfnString = argv[3];
    const char *errorsfn = errorsfnString.c_str();
    FILE *outfile = fopen (errorsfn, "w");
    buildPiecewiseLinear(reference, sa, outfile, printbuckets);
}
