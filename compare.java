/*
 * SAPLING: Suffix Array Piecewise Linear INdex for Genomics
 */
import java.util.*;
import java.io.*;
public class compare {
	static String refString;
	static boolean replace = true; // whether or not to replace all non-base characters with A's
	static int alpha = 2;
	static int k = 21;
	static int buckets = 15;
	static double mostThreshold = 0.95;
	static char[] reference;
	static int[] sa; //sa[i] is the location in the suffix array where character i in reference appears
	static int[] rev;
	static int n;
	static int[] errors;
	static SuffixArray lsa;
	static ArrayList<Integer> overs, unders;
	static int[] vals;
	static int maxOver, maxUnder, meanError;
	static int mostOver, mostUnder;
	static int[] llcp, rlcp;
public static void main(String[] args) throws IOException
{
	vals = new int[256];
	vals['A'] = (1<<alpha)-4;
	vals['C'] = (1<<alpha)-3;
	vals['G'] = (1<<alpha)-2;
	vals['T'] = (1<<alpha)-1;
	vals['N'] = (1<<alpha)-5;
	String fn = "/home/mkirsche/work/lis/chr22.fa";
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	StringBuilder sb = new StringBuilder("");
	input.nextLine();
	while(input.hasNext())
	{
		String cur = input.nextLine();
		if(!cur.startsWith(">"))
			sb.append(cur);
	}
	String s = sb.toString();
	if(replace) s = replace(s);
	refString = s;
	n = s.length();
	reference = s.toCharArray();
	System.out.println("Building suffix array of length " + s.length());
	lsa = new SuffixArray(s);
	sa = new int[n]; 
	rev = new int[n];
	for(int i = 0; i<n; i++) rev[sa[i] = lsa.inv[i]] = i;
	System.out.println("Built suffix array");
	buildPiecewiseLinear(s, lsa.inv);
	initializeLCPs();
	System.out.println("Mean error: " + meanError);
	System.out.println("Prediction range: (" + -maxUnder + ", " + maxOver + ")");
	System.out.println(mostThreshold + " of over-predictions are within " + mostOver);
	System.out.println(mostThreshold + " of under-predictions are within " + mostUnder);
	Random r = new Random(42);
	int numQueries = 5000000;
	char[][] queries = new char[numQueries][];
	long[] kmers = new long[numQueries];
	for(int i = 0; i<numQueries; i++)
	{
		int idx = r.nextInt(n - k);
		String query = s.substring(idx, idx + k);
		queries[i] = query.toCharArray();
		kmers[i] = kmerize(queries[i]);
	}
	// Run piece-wise linear test
	int[] plAnswers = new int[numQueries];
	long time = System.currentTimeMillis();
	for(int i = 0; i<numQueries; i++)
	{
		plAnswers[i] = plQuery(queries[i], kmers[i]);
	}
	long curTime = System.currentTimeMillis();
	long diff = curTime - time;
	System.out.println("Piecewise linear time: " + diff);
	
	// Check the answers for piecewise linear 
	int countCorrect = 0;
	for(int i = 0; i<numQueries; i++)
	{
		if(plAnswers[i] + k <= n && (new String(queries[i])).equals(s.substring(plAnswers[i], plAnswers[i]+k)))
		{
			countCorrect++;
		}
	}
	System.out.println("Piecewise linear correctness: " + countCorrect + " out of " + numQueries);
	
	//Run binary search test
	int[] bAnswers = new int[numQueries];
	time = System.currentTimeMillis();
	for(int i = 0; i<numQueries; i++)
	{
		//System.out.println(new String(queries[i]));
		bAnswers[i] = bQuery(queries[i]);
	}
	curTime = System.currentTimeMillis();
	diff = curTime - time;
	System.out.println("Binary search time: " + diff);
	
	// Check the answers for binary search 
	countCorrect = 0;
	for(int i = 0; i<numQueries; i++)
	{
		if(bAnswers[i] + k <= n && (new String(queries[i])).equals(s.substring(bAnswers[i], bAnswers[i]+k)))
		{
			countCorrect++;
		}
	}
	System.out.println("Binary search correctness: " + countCorrect + " out of " + numQueries);
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
	return rev[fancyBinarySearch(s, 0, n-k, loLcp, hiLcp)];
}

/*
 * Get the position in the reference of a query string using piecewise linear structure
 */
static int plQuery(char[] s)
{
	return plQuery(s, kmerize(s));
}

/*
 * Get the position in the reference of a query string using piecewise linear structure
 */
static int plQuery(char[] s, long kmer)
{
	int predicted = queryPiecewiseLinear(kmer); // Predicted position in suffix array
	int idx = rev[predicted]; // Actual string position where we predict it to be
	int lcp = getLcp(idx, s, 0);
	if(lcp == s.length) return idx;
	int lo, hi;
	int loLcp = -1, hiLcp = -1;
	if(lcp + idx == n || s[lcp] > reference[idx+lcp])
	{
		// Suffix is smaller then query - look farther right
		lo = predicted;
		hi = Math.min(n-1, predicted+mostOver); // Over-prediction which the actual position is highly likely to not exceed
		int hiIdx = rev[hi]; // String index corresponding to over-prediction
		int oLcp = getLcp(hiIdx, s, 0); // LCP between over-prediction suffix and query
		if(oLcp == s.length) return hiIdx; // Over-prediction happened to be exactly right
		if(oLcp + hiIdx == n || s[oLcp] > reference[hiIdx+oLcp])
		{
			// bad case: over-prediction still not high enough
			lo = hi;
			loLcp = oLcp;
			hi = Math.min(n-1, predicted + maxOver);
			hiIdx = rev[hi];
			oLcp = getLcp(hiIdx, s, 0);
			if(oLcp == s.length) return hiIdx;
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
		lo = Math.max(0, predicted-mostUnder);
		hi = predicted;
		int loIdx = rev[lo];
		int oLcp = getLcp(loIdx, s, 0); // LCP between under-prediction suffix and query
		if(oLcp == s.length) return loIdx; // Under-prediction happened to be exactly right
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
			lo = Math.max(0, predicted - maxUnder);
			loIdx = rev[lo];
			oLcp = getLcp(loIdx, s, 0);
			if(oLcp == s.length) return loIdx;
			loLcp = oLcp;
		}
	}
	return rev[binarySearch(s, lo, hi, loLcp, hiLcp)];
}
/*
 * Encode the first k characters of an an ACGT string into a 64-bit integer
 */
static long kmerize(char[] s)
{
	long kmer = 0;
	//if(alpha == 3 && ((k&1) == 0)) for(int i = 0; i<k; i+= 2) kmer = (kmer << 5) | (vals[s[i]] * 5 +vals[s[i+1]]);
	for(int i = 0; i<k; i++) kmer = (kmer << alpha) | vals[s[i]];
	return kmer;
}
/*
 * Replace all non-base characters in s with 'A'
 */
static String replace(String s)
{
	int n = s.length();
	char[] res = new char[n];
	for(int i = 0; i<n; i++)
	{
		char c = s.charAt(i);
		if(c == 'A' || c == 'C' || c == 'G' || c == 'T') res[i] = c;
		else res[i] = 'A';
	}
	return new String(res);
}
/*
 * Binary search.  We know:
 * Actual position in suffix array is in (lo, hi)
 * LCP between suffix corresponding to position lo in suffix array with query is loLcp
 * LCP between suffix corresponding to position hi in suffix array with query is hiLcp
 */
static int binarySearch(char[] s, int lo, int hi, int loLcp, int hiLcp)
{
	// Base case
	if(hi == lo + 2) return lo + 1;
	
	// TODO make this fancier and get rid of log factor
	int mid = (lo + hi) >> 1;
	int idx = rev[mid];
	int nLcp = getLcp(idx, s, Math.min(loLcp,  hiLcp));
	if(nLcp == s.length) return mid;
	if(nLcp + idx == n || s[nLcp] > reference[idx+nLcp])
	{
		// suffix too small - search right half
		return binarySearch(s, mid, hi, nLcp, hiLcp);
	}
	else
	{
		// suffix too big - search left half
		return binarySearch(s, lo, mid, loLcp, nLcp);
	}
}
/*
 * Fancy binary search.  We know:
 * Actual position in suffix array is in (lo, hi)
 * LCP between suffix corresponding to position lo in suffix array with query is loLcp
 * LCP between suffix corresponding to position hi in suffix array with query is hiLcp
 */
static int fancyBinarySearch(char[] s, int lo, int hi, int loLcp, int hiLcp)
{
	// Base case
	if(hi == lo+1) return -1;
	if(hi == lo + 2) return lo + 1;
	
	// TODO make this fancier and get rid of log factor
	int mid = (lo + hi) >> 1;
		
	if(lo<hi-1)
	{
		//System.out.println(lo+" "+hi+" "+refString.substring(rev[mid], Math.min(n, rev[mid]+20)));
		//System.out.println(refString.substring(rev[lo], Math.min(n, rev[lo]+20))+" "+refString.substring(rev[hi], Math.min(n, rev[hi]+20)));
		//System.out.println(loLcp+" "+hiLcp+" "+llcp[mid]+" "+rlcp[mid]);
	}
		
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
 * Fancy binary search modified to work when the initial call to the binary search can vary - currently too slow
 */
static int dynamicFancyBinarySearch(char[] s, int lo, int hi, int loLcp, int hiLcp)
{
	// Base case
	if(hi == lo + 2) return lo + 1;
		
	// TODO make this fancier and get rid of log factor
	int mid = (lo + hi) >> 1;
	if(loLcp >= hiLcp)
	{
		int llcp = lsa.lcp(rev[lo], rev[mid]);
		if(llcp > loLcp)
		{
			return dynamicFancyBinarySearch(s, mid, hi, loLcp, hiLcp);
		}
		else if(llcp == loLcp)
		{
			int idx = rev[mid];
			int nLcp = getLcp(idx, s, loLcp);
			if(nLcp == s.length) return mid;
			if(nLcp + idx == n || s[nLcp] > reference[idx+nLcp])
			{
				// suffix too small - search right half
				return dynamicFancyBinarySearch(s, mid, hi, nLcp, hiLcp);
			}
			else
			{
				// suffix too big - search left half
				return dynamicFancyBinarySearch(s, lo, mid, loLcp, nLcp);
			}
		}
		else
		{
			return dynamicFancyBinarySearch(s, lo, mid, loLcp, llcp);
		}
	}
	else
	{
		int rlcp = lsa.lcp(rev[hi], rev[mid]);
		if(rlcp > hiLcp)
		{
			return dynamicFancyBinarySearch(s, lo, mid, loLcp, hiLcp);
		}
		else if(rlcp == hiLcp)
		{
			int idx = rev[mid];
			int nLcp = getLcp(idx, s, hiLcp);
			if(nLcp == s.length) return mid;
			if(nLcp + idx == n || s[nLcp] > reference[idx+nLcp])
			{
				// suffix too small - search right half
				return dynamicFancyBinarySearch(s, mid, hi, nLcp, hiLcp);
			}
			else
			{
				// suffix too big - search left half
				return dynamicFancyBinarySearch(s, lo, mid, loLcp, nLcp);
			}
		}
		else
		{
			return dynamicFancyBinarySearch(s, mid, hi, rlcp, hiLcp);
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
/*
 * Code for the piecewise linear predictor of position in suffix array given a long representing a k-mer
 */
static long[] xlist, ylist;
static void buildPiecewiseLinear(String s, int[] sa)
{
	ArrayList<Long> xs = new ArrayList<Long>();
	ArrayList<Integer> ys = new ArrayList<Integer>();
	int[] vals = new int[256];
	for(int i = 0; i+k<=s.length(); i++)
	{
		long hash = kmerize(s.substring(i, i+k).toCharArray());
		xs.add(hash);
		ys.add(sa[i]);
	}
	xlist = new long[(1<<buckets)+1];
	ylist = new long[(1<<buckets)+1];
	Arrays.fill(xlist, -1);
	for(int i = 0; i<xs.size(); i++)
	{
		long x = xs.get(i), y = ys.get(i);
		int bucket = 0;
		bucket = (int)(x >> (alpha*k - buckets));
		if(xlist[bucket] == -1 || xlist[bucket] > x)
		{
			xlist[bucket] = x;
			ylist[bucket] = y;
		}
		if(x > xlist[xlist.length-1])
		{
			xlist[xlist.length-1] = x;
			ylist[ylist.length-1] = y;
		}
	}
	errors = new int[xs.size()];
	overs = new ArrayList<Integer>();
	unders = new ArrayList<Integer>();
	for(int i = 0; i<xs.size(); i++)
	{
		errors[i] = ys.get(i) - queryPiecewiseLinear(xs.get(i));
		if(errors[i] > 0) overs.add(errors[i]);
		else unders.add(-errors[i]);
	}
	errorStats();
}
static void errorStats()
{
	maxUnder = 0;
	maxOver = 0;
	long tot = 0;
	int n = errors.length;
	for(int i = 0; i<n; i++)
	{
		maxUnder = Math.max(-errors[i], maxUnder);
		maxOver = Math.max(errors[i], maxOver);
		tot += Math.abs(errors[i]);
	}
	meanError = (int)(.5 + tot / n);
	Collections.sort(overs);
	Collections.sort(unders);
	mostOver = overs.get((int)(mostThreshold * overs.size()));
	mostUnder = unders.get((int)(mostThreshold * unders.size()));
}
static int queryPiecewiseLinear(long x)
{
	int bucket = 0;
	//if(alpha == 3 && ((k&1) == 0)) bucket = (int)(x >> (5*(k>>1) - buckets));
	bucket = (int)(x >> (alpha*k - buckets));
	long xlo = xlist[bucket];
	long xhi = xlist[bucket+1];
	long ylo = ylist[bucket];
	long yhi = ylist[bucket+1];
	long predict = (long)(.5 + ylo + (yhi - ylo) * (x - xlo) * 1. / (xhi - xlo));
	return (int)predict;
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
