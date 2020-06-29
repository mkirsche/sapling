/*
 * An example of how to build and query Sapling.
 * This script also gives the timing of the Sapling queries in isolation, so is useful for benchmarking.
 */

#include <chrono>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include "sapling_api.h"

using namespace std::chrono;

int numBuckets = -1;
int maxMem = -1;
int k = -1;
int numQueries = 5000000;
string errorFnString = "";
string saFnString = "", saplingFnString = "";
int queryLength = -1;

Sapling sap;

void run_experiment(int queryLength);

int main(int argc, char **argv)
{
    if(argc < 2)
    {
        cout << "Usage: " << argv[0] << " <genome file> [saFn=<suffix array file>] [sapFn=<sapling file>] [nb=<log number of buckets>] [maxMem=<max number of buckets will be (genome size)/val>] [k=<k>] [nq=<number of queries>] [errFn=<errors file if outputting them>] [qLen=<query length>]" << endl;
        return 0;
    }

    string refFnString = argv[1];
    saFnString = refFnString + ".sa";
    saplingFnString = refFnString + ".sap";

    // Parse the command line arguments
    for(int i = 2; i<argc; i++)
    {
      string cur = argv[i];
      size_t eqPos = cur.find("=");
      if(eqPos != string::npos)
      {
        string arg = cur.substr(0, eqPos);
        string val = cur.substr(eqPos + 1);
        if(arg.compare("saFn") == 0)
        {
          saFnString = val;
        }
        if(arg.compare("sapFn") == 0)
        {
          saplingFnString = val;
        }
        if(arg.compare("errFn") == 0)
        {
          errorFnString = val;
        }
        if(arg.compare("nb") == 0)
        {
          numBuckets = stoi(val);
        }
        if(arg.compare("k") == 0)
        {
          k = stoi(val);
        }
        if(arg.compare("nq") == 0)
        {
          numQueries = stoi(val);
        }
        if(arg.compare("maxMem") == 0)
        {
          maxMem = stoi(val);
        }
		if(arg.compare("qLen") == 0)
        {
          queryLength = stoi(val);
        }
      }
    }

    // Build the Sapling data structure
    sap = Sapling(refFnString, saFnString, saplingFnString, numBuckets, maxMem, k, errorFnString);
    
    cout << "Testing Sapling" << endl;

	if(queryLength == -1)
	{
		run_experiment(sap.k-10);
		run_experiment(sap.k);
		run_experiment(sap.k + 10);
		run_experiment(sap.k + 20);
		run_experiment(sap.k + 30);
		run_experiment(sap.k + 80);
	}
    else
	{
		run_experiment(queryLength);
	}
}

void run_experiment(int queryLength)
{	
    cout << "Running experiment to search for " << queryLength << "-mers" << endl;
	// Create queries as random kmers from the genome
    vector<string> queries(numQueries, "");
    vector<long long> kmers = vector<long long>(numQueries, 0);
    vector<size_t> idxs = vector<size_t>(numQueries, 0);
    for(int i = 0; i<numQueries; i++)
    {
        idxs[i] = rand() % (sap.n - queryLength);
    	queries[i] = sap.reference.substr(idxs[i], queryLength);
    	kmers[i] = sap.kmerizeAdjusted(queryLength, queries[i]);
    }
	// Write the queries to a file
    FILE *outfile = fopen ("queries.out", "w");
    for(int i = 0; i < numQueries; i++)
    {
        fprintf(outfile, "@read%d\n", i+1);
        fprintf(outfile, "%s\n", queries[i].c_str());
        fprintf(outfile, "+\n");
        
        for(int j = 0; j<queryLength; j++) fprintf(outfile, "9");
        fprintf(outfile, "\n");
    }
    cout << "Constructed queries" << endl;
    
    // Run piece-wise linear test
    vector<long long> plAnswers(numQueries, 0);
    auto start = std::chrono::system_clock::now();
    for(int i = 0; i<numQueries; i++)
    {
        plAnswers[i] = sap.plQuery(queries[i].substr(0, queryLength), kmers[i], queries[i].length());
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "Piecewise linear time: " << elapsed_seconds.count() << endl;

    // Check the answers
    int countCorrect = 0;
    for(int i = 0; i<numQueries; i++)
    {
        if(plAnswers[i] == -1) continue;
        if(plAnswers[i] + (long long)queryLength <= (long long)sap.n && queries[i] == sap.reference.substr(plAnswers[i], queryLength))
        {
            countCorrect++;
        }
//else cout << idxs[i] << " " << plAnswers[i] << " " << queries[i] << " " << sap.reference.substr(plAnswers[i], queryLength) << endl;
    }
	cout << "Piecewise linear correctness: " << countCorrect << " out of " << numQueries << endl;
}
