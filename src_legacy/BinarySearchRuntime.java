/*
 * Measures the runtime of a binary search on a reference of different lengths
 */
import java.util.*;

import java.io.*;
public class BinarySearchRuntime {
	static int[] rev;
	static int alpha = 2;
	static int[] vals;
	static int[] llcp, rlcp;
	static char[] reference;
	static int n;
	static int k = 15;
	static SuffixArray lsa;
public static void main(String[] args) throws IOException
{
	vals = new int[256];
	vals['A'] = 0;
	vals['C'] = 1;
	vals['G'] = 2;
	vals['T'] = 3;
	String fn = "/home/mkirsche/work/lis/yeast.fa";
	PrintWriter out = new PrintWriter(new File("binarysearchruntime.txt"));
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	int maxLog = 23;
	int length = 0;
	char[] wholeRef = new char[1<<maxLog];
	while(input.hasNext() && length < (1<<maxLog))
	{
		String s = input.nextLine();
		if(s.startsWith(">")) continue;
		char[] cur = s.toCharArray();
		for(char c : cur)
		{
			if(length == 1<<(maxLog)) break;
			if(c == 'A' || c == 'C' || c == 'G' || c == 'T') wholeRef[length++] = c;
		}
	}
	String refString = new String(wholeRef);
	for(int log = 5; log <= maxLog; log++)
	{
		System.out.println(log);
		String cur = refString.substring(0,  1<<log);
		reference = cur.toCharArray();
		n = reference.length;
		lsa = new SuffixArray(cur);
		rev = new int[n];
		for(int i = 0; i<n; i++) rev[lsa.inv[i]] = i;
		initializeLCPs();
		Random r = new Random(42);
		double tot = 0;
		int iters = 20;
		for(int iter = 0; iter < iters; iter++)
		{
			int numQueries = 5000000;
			char[][] queries = new char[numQueries][];
			long[] kmers = new long[numQueries];
			//System.out.println(Arrays.toString(rev));
			for(int i = 0; i<numQueries; i++)
			{
				int idx = r.nextInt(n - k);
				//System.out.println(idx);
				String query = cur.substring(idx, idx + k);
				queries[i] = query.toCharArray();
				kmers[i] = kmerize(queries[i]);
			}
			int[] bAnswers = new int[numQueries];
			long time = System.currentTimeMillis();
			for(int i = 0; i<numQueries; i++)
			{
				bAnswers[i] = bQuery(queries[i]);
			}
			long curTime = System.currentTimeMillis();
			long diff = curTime - time;
			tot += diff;
		}
		out.println(log+" "+tot * 1.0 / iters);
		System.out.println("Runtime: " + tot * 1.0 / iters);
	}
	out.close();
}

static void initializeLCPs()
{
	llcp = new int[n];
	rlcp = new int[n];
	fillLCPs(0, n-k);
}
static void fillLCPs(int lo, int hi)
{
	if(hi - lo <= 2) return;
	int mid = (lo + hi) >> 1;
	llcp[mid] = lsa.lcp(rev[lo], rev[mid]);
	rlcp[mid] = lsa.lcp(rev[mid], rev[hi]);
	fillLCPs(lo, mid);
	fillLCPs(mid, hi);
}

/*
 * Get the position in the reference of a query string using piecewise binary search
 */
static int bQuery(char[] s)
{
	int loLcp = getLcp(rev[0], s, 0);
	if(loLcp == s.length) return rev[0];
	int hiLcp = getLcp(rev[n-k], s, 0);
	if(hiLcp == s.length) return rev[n-k];
	//System.out.println("Querying " + new String(s));
	int pos = fancyBinarySearch(s, 0, n-k, loLcp, hiLcp);
	//System.out.println(pos);
	return rev[pos];
}

/*
 * Encode the first k characters of an an ACGT string into a 64-bit integer
 */
static long kmerize(char[] s)
{
	long kmer = 0;
	for(int i = 0; i<k; i++) kmer = (kmer << alpha) | vals[s[i]];
	return kmer;
}

/*
 * Fancy binary search.  We know:
 * Actual position in suffix array is in (lo, hi)
 * LCP between suffix corresponding to position lo in suffix array with query is loLcp
 * LCP between suffix corresponding to position hi in suffix array with query is hiLcp
 */
static int fancyBinarySearch(char[] s, int lo, int hi, int loLcp, int hiLcp)
{
	//System.out.println(lo+" "+hi+" "+loLcp+" "+hiLcp);
	// Base case
	if(hi == lo+1) return -1;
	if(hi == lo + 2) return lo + 1;
	
	int mid = (lo + hi) >> 1;
		
	if(loLcp >= hiLcp)
	{
		if(llcp[mid] > loLcp)
		{
			return fancyBinarySearch(s, mid, hi, loLcp, hiLcp);
		}
		else if(llcp[mid] == loLcp)
		{
			int idx = rev[mid];
			int nLcp = getLcp(idx, s, loLcp);
			if(nLcp == s.length) return mid;
			if(nLcp + idx == n || s[nLcp] > reference[idx+nLcp])
			{
				// suffix too small - search right half
				return fancyBinarySearch(s, mid, hi, nLcp, hiLcp);
			}
			else
			{
				// suffix too big - search left half
				return fancyBinarySearch(s, lo, mid, loLcp, nLcp);
			}
		}
		else
		{
			return fancyBinarySearch(s, lo, mid, loLcp, llcp[mid]);
		}
	}
	else
	{
		if(rlcp[mid] > hiLcp)
		{
			return fancyBinarySearch(s, lo, mid, loLcp, hiLcp);
		}
		else if(rlcp[mid] == hiLcp)
		{
			int idx = rev[mid];
			int nLcp = getLcp(idx, s, hiLcp);
			if(nLcp == s.length) return mid;
			if(nLcp + idx == n || s[nLcp] > reference[idx+nLcp])
			{
				// suffix too small - search right half
				return fancyBinarySearch(s, mid, hi, nLcp, hiLcp);
			}
			else
			{
				// suffix too big - search left half
				return fancyBinarySearch(s, lo, mid, loLcp, nLcp);
			}
		}
		else
		{
			return fancyBinarySearch(s, mid, hi, rlcp[mid], hiLcp);
		}
	}
}

/*
 * Returns the lcp of the suffix of reference beginning at position idx and the string s
 */
static int getLcp(int idx, char[] s, int start)
{
	int i = start;
	for(; i<s.length && idx+i < n; i++) if(s[i] != reference[idx+i]) return i;
	return i; // whole thing matches
}

public static class SuffixArray {
    public RMQ rmq;
    public int[] str;
    public int[] idx;
    public int[] inv;
    public int[] lcp;

    public int length;
    public int letters;
    
    public static int[] toIntArray(String str) {
        int[] arr = new int[str.length() + 3];
        for (int i = 0; i < str.length(); i++) {
            arr[i] = str.charAt(i) - 'A' + 1;
        }
        return arr;
    }
    
    public SuffixArray(String dna)
    {
    	this(dna.replaceAll("C", "B").replaceAll("G", "C").replaceAll("T", "D").replaceAll("N", "E"), 1<<alpha);
    }

    public SuffixArray(String str, int letters) {
        // Assumes that string starts at 'a' and goes up to (char) 'a' + letters
        this(toIntArray(str), str.length(), letters + 1);
    }

    public SuffixArray(int[] str, int length, int letters) {
        this.length = length;
        this.letters = letters;
        this.str = str;
        this.idx = new int[length];

        build(str, idx, length, letters);

        this.inv = new int[length];
        for (int i = 0; i < length; i++) {
            inv[idx[i]] = i;
        }

        this.lcp = getLCP();
        this.rmq = new RMQ(this.lcp);
    }

    public void build(int[] str, int[] sa, int n, int letters) {
        int n0 = (n + 2) / 3, n2 = n / 3, n02 = n0 + n2, delta = n0 - (n + 1) / 3;

        int[] sa12 = new int[n02 + 3];
        int[] s12 = new int[n02 + 3];

        int[] sa0 = new int[n0];
        int[] s0 = new int[n0];

        // Sorting two thirds
        for (int i = 0, j = 0; i < n + delta; i++) {
            if (i % 3 > 0) {
                s12[j++] = i;
            }
        }

        radixPass(s12, sa12, str, 2, n02, letters);
        radixPass(sa12, s12, str, 1, n02, letters);
        radixPass(s12, sa12, str, 0, n02, letters);

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
            build(s12, sa12, n02, name);
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
        radixPass(s0, sa0, str, 0, n0, letters);

        // Merge
        for (int p = 0, t = delta, k = 0; k < n; k++) {
            int i = sa12[t] < n0 ? sa12[t] * 3 + 1 : (sa12[t] - n0) * 3 + 2;
            int j = sa0[p];
            if (sa12[t] < n0 ?
                    leq(str[i], s12[sa12[t] + n0], str[j], s12[j / 3]) :
                    leq(str[i], str[i + 1], s12[sa12[t] - n0 + 1],
                            str[j], str[j + 1], s12[j / 3 + n0])) {
                sa[k] = i;
                if (++t == n02) {
                    for (k++; p < n0; p++, k++) {
                        sa[k] = sa0[p];
                    }
                }
            } else {
                sa[k] = j;
                if (++p == n0) {
                    for (k++; t < n02; t++, k++) {
                        sa[k] = sa12[t] < n0 ? sa12[t] * 3 + 1 : (sa12[t] - n0) * 3 + 2;
                    }
                }
            }
        }
    }

    public boolean leq(int a1, int a2, int b1, int b2) {
        return (a1 < b1 || (a1 == b1 && a2 < b2));
    }
    public boolean leq(int a1, int a2, int a3, int b1, int b2, int b3) {
        return (a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3)));
    }

    public void radixPass(int[] a, int[] b, int[] ref, int offset, int n, int letters) {
        int[] cnt = new int[letters + 1];
        for (int i = 0; i < n; i++) {
            cnt[ref[a[i] + offset]]++;
        }
        for (int i = 0, sum = 0; i <= letters; i++) {
            int t = cnt[i];
            cnt[i] = sum;
            sum += t;
        }
        for (int i = 0; i < n; i++) {
            b[cnt[ref[a[i] + offset]]++] = a[i];
        }
    }

    public int[] getLCP() {
        int[] lcp = new int[length - 1];
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

    public int lcp(int aIdx, int bIdx) {
        if (aIdx == bIdx) {
            return length - aIdx;
        }
        int x = inv[aIdx];
        int y = inv[bIdx];
        return rmq.query(Math.min(x, y), Math.max(x, y) - 1);
    }
    static class RMQ
    {
    	RMQ(int[] a)
    	{
    		this.a = a;
    		preprocess();
    	}
    	int[][] rmq;
    	int[] a;
    	int[] memo;
    	int log(int x)
    	{
    		if(memo[x] != -1) return memo[x];
    		return memo[x] = Integer.numberOfTrailingZeros(Integer.highestOneBit(x));
    	}
    	void preprocess()
    	{
    		int n = a.length;
    		memo = new int[n+1];
    		Arrays.fill(memo, -1);
    		rmq = new int[n][log(n)+1];
    		for(int i = 0; i<n; i++) rmq[i][0] = i;
    		for(int j = 1; (1<<j) <= n; j++)
    			for(int i = 0; i + (1<<j) <= n; i++)
    				if(a[rmq[i][j-1]] < a[rmq[i+(1<<(j-1))][j-1]])
    					rmq[i][j] = rmq[i][j-1];
    				else rmq[i][j] = rmq[i+(1<<(j-1))][j-1];
    	}
    	int query(int i, int j)
    	{
    		int k = log(j - i + 1);
    		return Math.min(a[rmq[i][k]], a[rmq[j-(1<<k)+1][k]]);
    	}
    }
}

}
