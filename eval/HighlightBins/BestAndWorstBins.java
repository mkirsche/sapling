/*
 * Outputs points from some of the best and worst bins
 */

import java.util.*;
import java.io.*;
public class BestAndWorstBins {
	static int k = 21;
	static int bins = 18;
public static void main(String[] args) throws Exception
{
	String fn = "/home/mkirsche/git/SaplingEval/celegans.errors";
	String ofn = "celegansBinHighlights.txt";
	int highlightCount = 10;
	int minBinSize = 50;
	if(args.length == 3)
	{
		fn = args[0];
		ofn = args[1];
		highlightCount = Integer.parseInt(args[2]);
	}
	else
	{
		System.out.println("Usage: java BestAndWorstBins errorsFile numBins");
	}
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	bins = Integer.parseInt(input.nextLine().split(" ")[0]);
	HashMap<Long, Long> sums = new HashMap<Long, Long>();
	HashMap<Long, Integer> counts = new HashMap<Long, Integer>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		String[] tokens = line.split(" ");
		if(tokens.length == 4)
		{
			long x = Long.parseLong(tokens[0]);
			long error = Long.parseLong(tokens[3]);
			long binSize = 1 << (2 * k - bins);
			long binId = x / binSize;
			if(!sums.containsKey(binId))
			{
				sums.put(binId, Math.abs(error));
				counts.put(binId, 1);
			}
			else
			{
				sums.put(binId, sums.get(binId) + Math.abs(error));
				counts.put(binId, counts.get(binId) + 1);
			}
		}
	}
	
	ArrayList<Double> means = new ArrayList<Double>();
	
	for(long x : sums.keySet())
	{
		if(counts.get(x) < minBinSize)
		{
			continue;
		}
		means.add(1.0 * sums.get(x) / counts.get(x));
	}
	
	Collections.sort(means);
	
	double lowerBound = (means.size() < 2*highlightCount) ? -1 : means.get(highlightCount - 1);
	double upperBound = (means.size() < 2*highlightCount) ? -1 : means.get(means.size() - highlightCount);
		
	HashMap<Long, Integer> goodKeys = new HashMap<Long, Integer>();
	HashSet<Long> lowMean = new HashSet<Long>(); 
		
	for(long x : counts.keySet())
	{
		if(counts.get(x) < minBinSize)
		{
			continue;
		}
		double curMean = 1.0 * sums.get(x) / counts.get(x);
		if(curMean <= lowerBound + 1e-9 || curMean >= upperBound - 1e-9)
		{
			if(curMean <= lowerBound + 1e-9)
			{
				lowMean.add(x);
			}
			goodKeys.put(x, goodKeys.size());
		}
	}
	
	input = new Scanner(new FileInputStream(new File(fn)));
	PrintWriter out = new PrintWriter(new File(ofn));
	out.println(goodKeys.size());
	while(input.hasNext())
	{
		String line = input.nextLine();
		String[] tokens = line.split(" ");
		if(tokens.length == 4)
		{
			long x = Long.parseLong(tokens[0]);
			long error = Long.parseLong(tokens[3]);
			long binSize = 1 << (2 * k - bins);
			long binId = x / binSize;
			if(goodKeys.containsKey(binId))
			{
				out.println(binId+" "+goodKeys.get(binId)+" "+error+" "+(lowMean.contains(binId) ? "1" : "0"));
			}
		}
	}
	
	input.close();
	out.close();
}
static long max(ArrayList<Long> data)
{
	long res = data.get(0);
	for(int i = 1; i<data.size(); i++)
	{
		res = Math.max(res, data.get(i));
	}
	return res;
}

/*
 * Calculates the overall median error given the lists of errors in each bin
 */
static double median(HashMap<Long, ArrayList<Long>> errors)
{
  TreeMap<Long, Integer> counts = new TreeMap<Long, Integer>();
  long totalCount = 0;
  for(long key : errors.keySet())
  {
    ArrayList<Long> cur = errors.get(key);
    totalCount += cur.size();
    for(long error : cur)
    {
      counts.put(error, 1 + (counts.containsKey(error) ? counts.get(error) : 0));
    }
  }
  long totalFrequency = 0;
  long last = 0;
  for(long error : counts.keySet())
  {
    totalFrequency += counts.get(error);
    if(totalFrequency * 2 > totalCount)
    {
      if(counts.get(error) == 1 && (totalFrequency-1) * 2 == totalCount) return (error + last) / 2.0;
      return error * 1.0;
    }
    last = error;
  }
  return 0;
}

/*
 * Calculates the overall median error given the lists of errors in each bin
 */
static double percentile(HashMap<Long, ArrayList<Long>> errors, double prop)
{
  TreeMap<Long, Integer> counts = new TreeMap<Long, Integer>();
  long totalCount = 0;
  for(long key : errors.keySet())
  {
    ArrayList<Long> cur = errors.get(key);
    totalCount += cur.size();
    for(long error : cur)
    {
      counts.put(error, 1 + (counts.containsKey(error) ? counts.get(error) : 0));
    }
  }
  long totalFrequency = 0;
  long last = 0;
  for(long error : counts.keySet())
  {
    totalFrequency += counts.get(error);
    if(totalFrequency > totalCount * prop)
    {
      return error * 1.0;
    }
    last = error;
  }
  return 0;
}

static double mean(ArrayList<Long> data)
{
	double res = 0;
	for(long x : data) res += x;
	return res / data.size();
}
static double median(ArrayList<Long> data)
{
	Collections.sort(data);
	if(data.size()%2 == 0)
	{
		return (data.get(data.size() / 2) + data.get(data.size() / 2 - 1)) / 2.0;
	}
	return data.get(data.size() / 2);
}
static double maxDouble(ArrayList<Double> data)
{
	double res = data.get(0);
	for(int i = 1; i<data.size(); i++)
	{
		res = Math.max(res, data.get(i));
	}
	return res;
}
static double meanDouble(ArrayList<Double> data)
{
	double res = 0;
	for(double x : data) res += x;
	return res / data.size();
}
static double medianDouble(ArrayList<Double> data)
{
	Collections.sort(data);
	if(data.size()%2 == 0)
	{
		return (data.get(data.size() / 2) + data.get(data.size() / 2 - 1)) / 2.0;
	}
	return data.get(data.size() / 2);
}
}