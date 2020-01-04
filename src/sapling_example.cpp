/*
 * An example of how to build and query Sapling.
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

int main(int argc, char **argv)
{
    if(argc < 2)
    {
        cout << "Usage: " << argv[0] << " <genome file> [saFn=<suffix array file>] [sapFn=<sapling file>] [nb=<log number of buckets>] [maxMem=<max number of buckets will be (genome size)/val>] [k=<k>] [nq=<number of queries>] [errFn=<errors file if outputting them>]" << endl;
        return 0;
    }

    string refFnString = argv[1];
    saFnString = refFnString + ".sa";
    saplingFnString = refFnString + ".sap";

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
      }
    }

   /* if(argc >= 4)
    {
      saFnString = argv[2];
      saplingFnString = argv[3];
    }

    if(argc >= 5)
    {
      numBuckets = atoi(argv[4]);
    }*/

    Sapling sap(refFnString, saFnString, saplingFnString, numBuckets, maxMem, k, errorFnString);
    
    cout << "Testing Sapling" << endl;

    // Create queries as random kmers from the genome
    vector<string> queries(numQueries, "");
    vector<long long> kmers = vector<long long>(numQueries, 0);
    for(int i = 0; i<numQueries; i++)
    {
    	size_t idx = rand() % (sap.n - sap.k);
    	queries[i] = sap.reference.substr(idx, sap.k);
    	kmers[i] = sap.kmerize(queries[i]);
    }
    
    // Write the queries to a file
    FILE *outfile = fopen ("queries.out", "w");
    for(int i = 0; i < numQueries; i++)
    {
        fprintf(outfile, "@read%d\n", i+1);
        fprintf(outfile, "%s\n", queries[i].c_str());
        fprintf(outfile, "+\n");
        
        for(int j = 0; j<sap.k; j++) fprintf(outfile, "9");
        fprintf(outfile, "\n");
    }
    cout << "Constructed queries" << endl;
    
    // Run piece-wise linear test
    vector<long long> plAnswers(numQueries, 0);
    auto start = std::chrono::system_clock::now();
    for(int i = 0; i<numQueries; i++)
    {
        plAnswers[i] = sap.plQuery(queries[i].substr(0, sap.k), kmers[i], queries[i].length());
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "Piecewise linear time: " << elapsed_seconds.count() << endl;

    // Check the answers
    int countCorrect = 0;
    for(int i = 0; i<numQueries; i++)
    {
        if(plAnswers[i] == -1) continue;
        if(plAnswers[i] + (long long)sap.k <= (long long)sap.n && queries[i] == sap.reference.substr(plAnswers[i], sap.k))
        {
            countCorrect++;
        }
    }
    cout << "Piecewise linear correctness: " << countCorrect << " out of " << numQueries << endl;
}
