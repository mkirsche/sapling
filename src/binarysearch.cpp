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

vector<size_t> rev; // the inverse of sa sa[rev[i]] = i for all i
vector<unsigned int> llcp, rlcp;

string reference;

int k = 21;
int alpha = 2;
size_t n;

SuffixArray lsa;

size_t getLcp(size_t idx, string s, size_t start, size_t length)
{
	size_t i = start;
	for(; i<length && idx+i < n; i++) if(s[i] != reference[idx+i]) return i;
	return i; // whole thing matches
}

/*
 * Binary search.  We know:
 * Actual position in suffix array is in (lo, hi)
 * LCP between suffix corresponding to position lo in suffix array with query is loLcp
 * LCP between suffix corresponding to position hi in suffix array with query is hiLcp
 */
size_t binarySearch(string &s, size_t lo, size_t hi, size_t loLcp, size_t hiLcp, int length)
{
  //cout << lo << " " << hi << " " << loLcp << " " << hiLcp << " " << length << endl;
	// Base case
	if(hi == lo + 2) return lo + 1;
	
	size_t mid = (lo + hi) >> 1;
	size_t idx = rev[mid];
	int nLcp = getLcp(idx, s, min(loLcp,  hiLcp), length);
	if(nLcp == length) return mid;
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

size_t calcRLCP(size_t lo, size_t hi);

size_t calcLLCP(size_t lo, size_t hi)
{
  if(hi == lo + 2) return lsa.lcp[lo];
  if(hi == lo + 1) return n;
	size_t mid = (lo + hi) >> 1;
	size_t res = min(calcLLCP(lo, mid), calcRLCP(lo, mid));
	llcp[mid] = (unsigned int)res;
	return res;
}

size_t calcRLCP(size_t lo, size_t hi)
{
  if(hi == lo + 2) return lsa.lcp[lo+1];
  if(hi == lo + 1) return lsa.lcp[lo];
  size_t mid = (lo + hi) >> 1;
	size_t res = min(calcLLCP(mid, hi), calcRLCP(mid, hi));
	rlcp[mid] = (unsigned int)res;
	return res;
}

static void initializeLCPs()
{
	llcp.resize(n);
	rlcp.resize(n);
	calcLLCP(0, n-k);
	calcRLCP(0, n-k);
}

size_t fancyBinarySearch(string s, size_t lo, size_t hi, size_t loLcp, size_t hiLcp, int length)
{
	// Base case
	if(hi == lo+1) return n+1;
	if(hi == lo + 2) return lo + 1;
	
	// TODO make this fancier and get rid of log factor
	size_t mid = (lo + hi) >> 1;
		
	if(loLcp >= hiLcp)
	{
		if(llcp[mid] > loLcp)
		{
			return fancyBinarySearch(s, mid, hi, loLcp, hiLcp, length);
		}
		else if(llcp[mid] == loLcp)
		{
			size_t idx = rev[mid];
			int nLcp = getLcp(idx, s, loLcp, length);
			if(nLcp == length) return mid;
			if(nLcp + idx == n || s[nLcp] > reference[idx+nLcp])
			{
				// suffix too small - search right half
				return fancyBinarySearch(s, mid, hi, nLcp, hiLcp, length);
			}
			else
			{
				// suffix too big - search left half
				return fancyBinarySearch(s, lo, mid, loLcp, nLcp, length);
			}
		}
		else
		{
			return fancyBinarySearch(s, lo, mid, loLcp, llcp[mid], length);
		}
	}
	else
	{
		if(rlcp[mid] > hiLcp)
		{
			return fancyBinarySearch(s, lo, mid, loLcp, hiLcp, length);
		}
		else if(rlcp[mid] == hiLcp)
		{
			size_t idx = rev[mid];
			int nLcp = getLcp(idx, s, hiLcp, length);
			if(nLcp == length) return mid;
			if(nLcp + idx == n || s[nLcp] > reference[idx+nLcp])
			{
				// suffix too small - search right half
				return fancyBinarySearch(s, mid, hi, nLcp, hiLcp, length);
			}
			else
			{
				// suffix too big - search left half
				return fancyBinarySearch(s, lo, mid, loLcp, nLcp, length);
			}
		}
		else
		{
			return fancyBinarySearch(s, mid, hi, rlcp[mid], hiLcp, length);
		}
	}
}

/*
 * Get the position in the reference of a query string using piecewise binary search
 */
size_t bQuery(string s)
{
	size_t loLcp = getLcp(rev[0], s, 0, s.length());
	if(loLcp == s.length()) return rev[0];
	size_t hiLcp = getLcp(rev[n-1], s, 0, s.length());
	if(hiLcp == s.length()) return rev[n-1];
	return rev[binarySearch(s, 0, n-1, loLcp, hiLcp, s.length())];
}

int main(int argc, char **argv)
{
    if(argc != 3)
    {
        cout << "Usage: " << argv[0] << " <genome> " << " <suffix array file> " << endl;
        return 0;
    }
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
    
    string fnString = argv[2];
    const char *fn = fnString.c_str();
    ifstream f(fn);
    if(f.good())
    {
        cout << "Reading suffix array from file" << endl;
        FILE *infile = fopen (fn, "rb");
        size_t size;
        fread(&size, sizeof(size_t), 1, infile);
        lsa = SuffixArray();
        lsa.inv.resize(size);
        fread(&lsa.inv[0], sizeof(size_t), size, infile);
        
        fread(&size, sizeof(size_t), 1, infile);
        lsa.lcp.resize(size);
        fread(&lsa.lcp[0], sizeof(size_t), size, infile);
        
        cout << "Loaded suffix array of size " << lsa.inv.size() << endl;
    }
    else
    {
        cout << "Building suffix array" << endl;
        lsa = sa_init3(reference, alpha);
        cout << "Writing suffix array to file" << endl;
        FILE *outfile = fopen (fn, "wb");
        size_t size = lsa.inv.size();
        fwrite(&size, sizeof(size_t), 1, outfile);
        fwrite(&lsa.inv[0], sizeof(size_t), size, outfile);
        size = lsa.lcp.size();
        fwrite(&size, sizeof(size_t), 1, outfile);
        fwrite(&lsa.lcp[0], sizeof(size_t), size, outfile);
        cout << "Built suffix array of size " << lsa.inv.size() << endl;
    }
    
    cout << "Initializing rev and sa" << endl;
    rev = vector<size_t>(n, 0);
    cout << "Filling rev and sa" << endl;
    for(size_t i = 0; i<n; i++) rev[lsa.inv[i]] = i;
    
    cout << "Initializing LCP arrays" << endl;
    initializeLCPs();
    
    cout << "Testing binary search" << endl;
    int numQueries = 5000000;
    vector<string> queries(numQueries, "");
    vector<size_t> indexes(numQueries, 0);
    for(int i = 0; i<numQueries; i++)
    {
    	size_t idx = rand() % (n - k);
      indexes[i] = idx;
    	queries[i] = reference.substr(idx, k);
    }
    cout << "Constructed queries" << endl;
    // Run binary search test
    vector<size_t> bAnswers(numQueries, 0);
    auto start = std::chrono::system_clock::now();
    for(int i = 0; i<numQueries; i++)
    {
        bAnswers[i] = bQuery(queries[i]);
    }
    
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "Binary search time: " << elapsed_seconds.count() << endl;
    
    // Check the answers
    int countCorrect = 0;
    for(int i = 0; i<numQueries; i++)
    {
        if(bAnswers[i] + k <= n && queries[i] == reference.substr(bAnswers[i], k))
        {
            countCorrect++;
        }
        else
        {
            cout << queries[i] << " " << indexes[i] << " " << lsa.inv[indexes[i]] << " " << bAnswers[i] << endl;
        }
    }
    cout <<"Binary search correctness: " << countCorrect << " out of " << numQueries << endl;
    return 0;
}
