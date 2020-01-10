import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Random;
public class SuffixArraySimulatedSequences {
	static Random rand;
public static void main(String[] args) throws Exception
{
	rand = new Random(1212121);
	int length = 100000;
	int k = 21;
	String[] tests = new String[6];
	tests = new String[] {
			simulate(length, new double[] {1, 1, 1, 1}),
			simulate(length, new double[] {2, 3, 3, 2}),
			simulate(length, new double[] {1, 3, 3, 1}),
			simulate(length, new double[] {1, 9, 9, 1}),
			simulate(length, new double[] {0, 1, 1, 0}),
			simulate(length, new double[] {0, 1, 0, 0}),
			simulateRepeat(length, "CT"),
			simulateRepeat(length, "CAT"),
			simulateRepeat(length, "ACGT"),
			simulateRepeat(length, "ACTTCA"),
			simulateRepeat(length, "ACGTCGTAGTACTACG"),
			simulateRepeat(length, "ATGCTAACTAGGTCGATCGTATGCTAACTTTGCTAAGCTTTCCCTACTGC"),
	};
	for(int testNum = 0; testNum < tests.length; testNum++)
	{
		String test = tests[testNum];
		int n = test.length();
		String ofn = "suffixArraySim" + (testNum + 1) + ".txt";
		PrintWriter out = new PrintWriter(new File(ofn));
		suffixArray(test);
		int[] sa = rank[last];
		int[] inv = new int[n];
		for(int i = 0; i<n; i++) inv[sa[i]] = i;
		for(int i = 0; i+k <= n; i++)
		{
			String substr = test.substring(i, i+k);
			int saPos = sa[i];
			out.println(i + " " + kmerize(substr) + " " + saPos);
		}
		out.close();
	}
}

/*
 * Gets the numeric value of a single base paiir (A = 0, C = 1, G = 2, T = 3)
 */
static int bpToInt(char c)
{
	if(c >= 'a' && c <= 'z') c += 'A' - 'a';
	if(c == 'A') return 0;
	if(c == 'C') return 1;
	if(c == 'G') return 2;
	return 3;
}

/*
 * Converts a k-mer string to an integer with 2k bits
 */
static long kmerize(String seq)
{
	long res = 0;
	for(int i = 0; i<seq.length(); i++)
	{
		char c = seq.charAt(i);
		res *= 4;
		res |= bpToInt(c);
	}
	return res;
}

/*
 * Simulates a random sequence with given length and nucleotide frequencies
 */
static char[] bps = new char[] {'A', 'C', 'G', 'T'};
static String simulate(int length, double[] fs)
{
	// Normalize nucleotide frequencies
	double sum = 0;
	for(double f : fs) sum += f;
	for(int i = 0; i<fs.length; i++) fs[i] /= sum;
	
	int n = length;
	char[] res = new char[n];
	for(int i = 0; i<n; i++)
	{
		double val = rand.nextDouble();
		res[i] = bps[bps.length-1];
		for(int j = 0; j<bps.length-1; j++)
		{
			val -= fs[j];
			if(val <= 0)
			{
				res[i] = bps[j];
				break;
			}
		}
	}
	return new String(res);
}

/*
 * Makes a sequence of a given length consisting of a given repeat
 */
static String simulateRepeat(int length, String rep)
{
	int n = length;
	char[] res = new char[n];
	for(int i = 0; i<n; i++) res[i] = rep.charAt(i%rep.length());
	return new String(res);
}

/*
 * Suffix Arrays - O(n*lg^2(n))
 * Computes array s.t. rank[i] has a ranking after the first 2^i characters are looked.
 * rank[last] has a ranking of all the substring
 * Example: banana -> rank[last] is {3,2,5,1,4,0} because "banana" is the 4th in sorted array, "anana" is third, and so on.
 */
static int[][] rank;
static int last;
static void suffixArray(String s)
{
	int n = s.length();
	rank = new int[(int)Math.ceil((Math.log(n)/Math.log(2)))+2][n];
	for(int i = 0; i<n; i++) rank[0][i] = s.charAt(i) - 'A';
	for(int S = 1, L = 1; L/2 < n; L *= 2, S++)
	{
		Suffix[] suff = new Suffix[n];
		for(int i = 0; i<n; i++) suff[i] = new Suffix(rank[S-1][i], i+L >= n ? -1 : rank[S-1][i+L], i);
		Arrays.sort(suff);
		for(int i = 1; i<n; i++) rank[S][suff[i].i] = (suff[i].compareTo(suff[i-1]) == 0) ? rank[S][suff[i-1].i] : i;
		last = S;
	}
}
//O(lg(n)) -> gives length of longest common prefix given two indices of suffixes
static int lcp(int x, int y)
{
	int n = rank[0].length, ret = 0;
	if(x == y) return n - x;
	for(int i = last; i>= 0; i--)
	{
		if(rank[i][x] == rank[i][y])
		{
			ret += (1<<i);
			x += (1<<i);
			y += (1<<i);
			if(x >= n || y >= n) return ret;
		}
	}
	return ret;
}
static class Suffix implements Comparable<Suffix>
{
	int pr, cr, i;
	Suffix(int p, int c, int ii)
	{
		pr = p; cr = c; i = ii;
	}
	public int compareTo(Suffix o) 
	{
		if(this.pr == o.pr) return this.cr - o.cr;
		return this.pr - o.pr;
	}
}
}
