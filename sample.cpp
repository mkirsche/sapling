#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include "sa.h"
#include "util.h"

int k = 21;
int alpha = 2;
int vals[256];

long long kmerize(const string& s)
{
    long long kmer = 0;
	for(int i = 0; i<k; i++) kmer = (kmer << alpha) | vals[s[i]];
	return kmer;
}

int main(int argc, char **argv)
{
    for(int i = 0; i<256; i++) vals[i] = 0;
    vals['A'] = (1<<alpha)-4;
    vals['C'] = (1<<alpha)-3;
    vals['G'] = (1<<alpha)-2;
    vals['T'] = (1<<alpha)-1;
    vals['N'] = (1<<alpha)-5;
    ifstream input(argv[1]);
    string cur;
    std::ostringstream out("");
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
    string reference = out.str();
    int n = reference.length();
    
    string fnString = argv[2];
    const char *fn = fnString.c_str();
    FILE *infile = fopen (fn, "rb");
    size_t size;
    fread(&size, sizeof(size_t), 1, infile);
    size_t x;
    for(size_t i = 0; i < size; i++)
    {
        fread(&x, sizeof(size_t), 1, infile);
        if((i%1000 == 999) && (i + k <= n))
        {
            cout << kmerize(reference.substr(i, k)) << ' ' << x << endl;
        }
    }
    return 0;
}
