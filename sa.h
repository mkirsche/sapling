/*
 * Suffix Array Implementation
 */
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
 
using namespace std;

struct RMQ
{
    vector<int> a;
    vector<vector<int> > rmq;
    RMQ(vector<int> aa)
    {
        a = aa;
        int n = aa.size();
        rmq.resize(n);
        for (int i = 0; i < n; ++i)
            rmq[i].resize(log(n)+1);
        for(int i = 0; i<n; i++) rmq[i][0] = i;
    		for(int j = 1; (1<<j) <= n; j++)
    			for(int i = 0; i + (1<<j) <= n; i++)
    				if(a[rmq[i][j-1]] < a[rmq[i+(1<<(j-1))][j-1]])
    					rmq[i][j] = rmq[i][j-1];
    				else rmq[i][j] = rmq[i+(1<<(j-1))][j-1];
    }
    RMQ(){}
    int log(int n)
    {
        int res = 0;
        while(n >= (1<<res)) res++;
        return res-1;
    }
    int query(int i, int j)
    {
        int k = log(j - i + 1);
    	return min(a[rmq[i][k]], a[rmq[j-(1<<k)+1][k]]);
    }
};

vector<int> toIntVector(string str) {
    vector<int> res(str.length() + 3, 0);
    for (int i = 0; i < str.length(); i++) {
        res[i] = str[i] - 'A' + 1;
    }
    res[str.length()] = res[str.length()+1] = res[str.length()+2] = 0;
    return res;
}
 
 struct SuffixArray {
    vector<int> str;
    vector<int> idx;
    vector<int> inv;
    vector<int> lcp;
    
    RMQ rmq;

    int length;
    int letters;

    void radixPass(vector<int> a, vector<int>* b, vector<int> ref, int offset, int n, int letters) {
        vector<int> cnt(letters+1, 0);
        for(int i = 0; i<letters+1; i++) cnt[i] = 0;
        for (int i = 0; i < n; i++) {
            cnt[ref[a[i] + offset]]++;
        }
        for (int i = 0, sum = 0; i <= letters; i++) {
            int t = cnt[i];
            cnt[i] = sum;
            sum += t;
        }
        for (int i = 0; i < n; i++) {
            (*b)[cnt[ref[a[i] + offset]]++] = a[i];
        }
    }

    void build(vector<int> str, vector<int>* sap, int n, int letters) {
        int n0 = (n + 2) / 3, n2 = n / 3, n02 = n0 + n2, delta = n0 - (n + 1) / 3;
        
        vector<int> sa12(n02 + 3, 0);
        vector<int> s12(n02 + 3, 0);
        
        vector<int> sa0(n0, 0);
        vector<int> s0(n0, 0);

        // Sorting two thirds
        for (int i = 0, j = 0; i < n + delta; i++) {
            if (i % 3 > 0) {
                s12[j++] = i;
            }
        }
        radixPass(s12, &sa12, str, 2, n02, letters);
        radixPass(sa12, &s12, str, 1, n02, letters);
        radixPass(s12, &sa12, str, 0, n02, letters);

        // Checking if the suffixes are sufficiently sorted
        int name = 0, c0 = -1, c1 = -1, c2 = -1;
        for (int i = 0; i < n02; i++) {
            if (str[sa12[i]] != c0 || str[sa12[i] + 1] != c1 || str[sa12[i] + 2] != c2) {
                name++;
                c0 = str[sa12[i]];
                c1 = str[sa12[i] + 1];
                c2 = str[sa12[i] + 2];
            }
            if (sa12[i] % 3 == 1) {
                s12[sa12[i] / 3] = name;
            } else {
                s12[sa12[i] / 3 + n0] = name;
            }
        }

        // Recursively sort if not, generate array if it is
        if (name < n02) {
            build(s12, &sa12, n02, name);
            for (int i = 0; i < n02; i++) {
                s12[sa12[i]] = i + 1;
            }
        } else {
            for (int i = 0; i < n02; i++) {
                sa12[s12[i] - 1] = i;
            }
        }

        // Sorting lone third
        for (int i = 0, j = 0; i < n02; i++) {
            if (sa12[i] < n0) {
                s0[j++] = 3 * sa12[i];
            }
        }
        radixPass(s0, &sa0, str, 0, n0, letters);

        // Merge
        for (int p = 0, t = delta, k = 0; k < n; k++) {
            int i = sa12[t] < n0 ? sa12[t] * 3 + 1 : (sa12[t] - n0) * 3 + 2;
            int j = sa0[p];
            if (sa12[t] < n0 ?
                    leq(str[i], s12[sa12[t] + n0], str[j], s12[j / 3]) :
                    leq(str[i], str[i + 1], s12[sa12[t] - n0 + 1],
                            str[j], str[j + 1], s12[j / 3 + n0])) {
                (*sap)[k] = i;
                if (++t == n02) {
                    for (k++; p < n0; p++, k++) {
                        (*sap)[k] = sa0[p];
                    }
                }
            } else {
                (*sap)[k] = j;
                if (++p == n0) {
                    for (k++; t < n02; t++, k++) {
                        (*sap)[k] = sa12[t] < n0 ? sa12[t] * 3 + 1 : (sa12[t] - n0) * 3 + 2;
                    }
                }
            }
        }
    }

    int leq(int a1, int a2, int b1, int b2) {
        return (a1 < b1 || (a1 == b1 && a2 < b2));
    }
    int leq(int a1, int a2, int a3, int b1, int b2, int b3) {
        return (a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3)));
    }

    vector<int> getLCP() {
        vector<int> lcp(length - 1, 0);
        int curr = 0;
        for (int i = 0; i < length; i++) {
            int k = inv[i];
            if (k < length - 1) {
                int j = idx[k + 1];
                while (i + curr < length && j + curr < length &&
                        str[i + curr] == str[j + curr]) {
                    curr++;
                }
                lcp[k] = curr;
            }
            if (curr > 0) {
                curr--;
            }
        }
        return lcp;
    }
    
    int queryLcp(int a, int b)
    {
        if(a == b) return length - a;
        int x = inv[a], y = inv[b];
        return rmq.query(min(x, y), max(x, y)-1);
    }
    
    SuffixArray(vector<int> st, int ln, int ltrs)
    {
        length = ln;
        letters = ltrs;
        str = st;
        idx.resize(length);
        build(str, &idx, length, letters);
        inv.resize(length);
        for (int i = 0; i < length; i++) {
            inv[idx[i]] = i;
        }

        lcp = getLCP();
        rmq = RMQ(lcp);
    }
    
    SuffixArray(){}
};

SuffixArray sa_init2(string str, int letters)
{
    return SuffixArray(toIntVector(str), str.length(), letters + 1);
}

SuffixArray sa_init3(string dna, int alpha)
{
    string rep = dna;
    std::replace(rep.begin(), rep.end(), 'C', 'B');
    std::replace(rep.begin(), rep.end(), 'G', 'C');
    std::replace(rep.begin(), rep.end(), 'T', 'D');
    std::replace(rep.begin(), rep.end(), 'N', 'E');
    return sa_init2(rep, 1<<alpha);
}

