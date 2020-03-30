/*
 * Suffix Array Implementation
 */
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#define myint unsigned int
 
using namespace std;

vector<size_t> toIntVector(string str) {
    vector<size_t> res(str.length() + 3, 0);
    for (size_t i = 0; i < str.length(); i++) {
        res[i] = (size_t)(str[i] - 'A' + 1);
    }
    res[str.length()] = res[str.length()+1] = res[str.length()+2] = (size_t)0;
    return res;
}
 
struct SuffixArray {

    vector<size_t> str;
    vector<size_t> idx;
    vector<size_t> inv;
    vector<size_t> lcp;
    
    size_t krmqk = 0;
    vector<myint> krmqb;

    // Say for each position in the lcp array, how far to the right the next lcp value less then k is which occurs at or after that position
    void krmq_init(size_t kk)
    {
      krmqk = kk;
      size_t n = lcp.size();
      krmqb.resize(n+1);
      krmqb[n] = 0;
      for(size_t i = n; i-->0 ;)
      {
          krmqb[i] = (lcp[i] < krmqk) ? 0 : (1 + krmqb[i+1]);
      }
    }

    // Whether suffixes indexed i and j in sorted order have LCP >= k
    // Assumes i < j
    int krmq_query(size_t i, size_t j)
    {
        //cout << i << " " << j << " " << krmqb[i] << endl;
        return (i > j) || ((size_t)krmqb[i] +  i > j);
    }

    // Whether suffixes with indices a and b in sorted order have lcp >= k
    int queryLcpK(size_t a, size_t b)
    {
        return krmq_query(min(a, b), max(a, b)-1);
    }

    size_t length = 0;
    size_t letters = 0;

    void radixPass(vector<size_t> a, vector<size_t>* b, vector<size_t> ref, size_t offset, size_t n, size_t letters)
    {
	    vector<size_t> cnt(letters+1, 0);
      for(size_t i = 0; i<letters+1; i++) cnt[i] = 0;
      for (size_t i = 0; i < n; i++)
      {
        cnt[(size_t)(ref[a[i] + offset])]++;
      }
      for (size_t i = 0, sum = 0; i <= letters; i++)
      {
        size_t t = cnt[i];
        cnt[i] = sum;
        sum += t;
      }
      for (size_t i = 0; i < n; i++)
      {
        (*b)[cnt[(size_t)(ref[a[i] + offset])]++] = a[i];
      }
    }

    void build(vector<size_t> str, vector<size_t>* sap, size_t n, size_t letters) 
    {
      size_t n0 = (n + 2) / 3, n2 = n / 3, n02 = n0 + n2, delta = n0 - (n + 1) / 3;
        
      vector<size_t> sa12(n02 + 3, 0);
      vector<size_t> s12(n02 + 3, 0);
        
      vector<size_t> sa0(n0, 0);
      vector<size_t> s0(n0, 0);

      // Sorting two thirds
      for (size_t i = 0, j = 0; i < n + delta; i++)
      {
          if (i % 3 > 0) 
          {
            s12[j++] = i;
          }
      }
      radixPass(s12, &sa12, str, 2, n02, letters);
      radixPass(sa12, &s12, str, 1, n02, letters);
      radixPass(s12, &sa12, str, 0, n02, letters);

      // Checking if the suffixes are sufficiently sorted
      size_t name = 0, c0 = -1, c1 = -1, c2 = -1;
      for (size_t i = 0; i < n02; i++)
      {
        if (str[sa12[i]] != c0 || str[sa12[i] + 1] != c1 || str[sa12[i] + 2] != c2)
        {
          name++;
          c0 = str[sa12[i]];
          c1 = str[sa12[i] + 1];
          c2 = str[sa12[i] + 2];
        }
        if (sa12[i] % 3 == 1)
        {
          s12[sa12[i] / 3] = name;
        }
        else
        {
          s12[sa12[i] / 3 + n0] = name;
        }
      }

      // Recursively sort if not, generate array if it is
      if (name < n02)
      {
        build(s12, &sa12, n02, name);
        for (size_t i = 0; i < n02; i++)
        {
          s12[sa12[i]] = i + 1;
        }
      } 
      else
      {
        for (size_t i = 0; i < n02; i++)
        {
          sa12[s12[i] - 1] = i;
        }
      }

      // Sorting lone third
      for (size_t i = 0, j = 0; i < n02; i++)
      {
        if (sa12[i] < n0)
        {
          s0[j++] = 3 * sa12[i];
        }
      }
      radixPass(s0, &sa0, str, 0, n0, letters);

      // Merge
      for (size_t p = 0, t = delta, k = 0; k < n; k++) 
      {
        size_t i = sa12[t] < n0 ? sa12[t] * 3 + 1 : (sa12[t] - n0) * 3 + 2;
        size_t j = sa0[p];
        if (sa12[t] < n0 ?
          leq(str[i], s12[sa12[t] + n0], str[j], s12[j / 3]) :
          leq(str[i], str[i + 1], s12[sa12[t] - n0 + 1],
          str[j], str[j + 1], s12[j / 3 + n0]))
        {
          (*sap)[k] = i;
          if (++t == n02) 
          {
            for (k++; p < n0; p++, k++)
            {
              (*sap)[k] = sa0[p];
            }
          }
        }
        else
        {
          (*sap)[k] = j;
          if (++p == n0)
          {
            for (k++; t < n02; t++, k++)
            {
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

    vector<size_t> getLCP() {
        vector<size_t> lcp(length - 1, 0);
        size_t curr = 0;
        for (size_t i = 0; i < length; i++) {
            size_t k = inv[i];
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

    int queryLcpFromSAPos(size_t a, size_t b)
	  {
		  return krmq_query(a, b-1);
	  }
    
    SuffixArray(vector<size_t> st, size_t ln, size_t ltrs)
    {
        length = ln;
        letters = ltrs;
        str = st;
        idx.resize(length);
        build(str, &idx, length, letters);
        inv.resize(length);
        for (size_t i = 0; i < length; i++) {
            inv[idx[i]] = i;
        }
        lcp = getLCP();
        vector<size_t>().swap(str);
        vector<size_t>().swap(idx);
    }
    
    SuffixArray(){}
};

SuffixArray sa_init2(string str, size_t letters)
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

