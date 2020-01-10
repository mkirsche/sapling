import java.util.*;
import java.io.*;
public class PerBinErrors {
	static int k = 21;
	static int bins = 18;
public static void main(String[] args) throws Exception
{
	String fn = "/home/mkirsche/git/SaplingEval/allerrors.txt";
	String ofn = "stats.txt";
	if(args.length == 2)
	{
		fn = args[0];
		ofn = args[1];
	}
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
  bins = Integer.parseInt(input.nextLine().split(" ")[0]);
	HashMap<Long, ArrayList<Long>> errors = new HashMap<Long, ArrayList<Long>>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		String[] tokens = line.split(" ");
		if(tokens.length == 4)
		{
			long x = Long.parseLong(tokens[0]);
			long error = Long.parseLong(tokens[3]);
			long binSize = 1L << (2 * k - bins);
			long binId = x / binSize;
			if(!errors.containsKey(binId))
			{
				errors.put(binId, new ArrayList<Long>());
			}
			errors.get(binId).add(Math.abs(error));
		}
	}
	
	PrintWriter out = new PrintWriter(new File(ofn));
	
	ArrayList<Long> maxes = new ArrayList<Long>();
	ArrayList<Double> means = new ArrayList<Double>();
	ArrayList<Double> medians = new ArrayList<Double>();
	
	for(long x : errors.keySet())
	{
		ArrayList<Long> cur = errors.get(x);
		maxes.add(max(cur));
		means.add(mean(cur));
		medians.add(median(cur));
		out.println(maxes.get(maxes.size()-1)+" "+medians.get(medians.size()-1)+" "+means.get(means.size()-1));
	}

  System.out.println("Median: " + median(errors));
  System.out.println("95th percentile: " + percentile(errors, .95));
	
	System.out.println("Max of maxes: " + max(maxes));
	System.out.println("Median of maxes: " + median(maxes));
	System.out.println("Mean of maxes: " + mean(maxes));
	Collections.sort(maxes);
	System.out.print(Math.min(20, maxes.size()) + " highest maxes: ");
	for(int i = Math.max(0, maxes.size() - 20); i < maxes.size(); i++)
	{
		System.out.print(maxes.get(i)  + " ");
	}
	System.out.println();
	System.out.println("Max of medians: " + maxDouble(medians));
	System.out.println("Median of medians: " + medianDouble(medians));
	System.out.println("Mean of medians: " + meanDouble(medians));
	Collections.sort(medians);
	System.out.print(Math.min(20, medians.size()) + " highest medians: ");
	for(int i = Math.max(0, medians.size() - 20); i < medians.size(); i++)
	{
		System.out.print(medians.get(i)  + " ");
	}
	System.out.println();
	System.out.println("Max of means: " + maxDouble(means));
	System.out.println("Median of means: " + medianDouble(means));
	System.out.println("Mean of means: " + meanDouble(means));
	Collections.sort(means);
	System.out.print(Math.min(20, means.size()) + " highest means: ");
	for(int i = Math.max(0, medians.size() - 20); i < means.size(); i++)
	{
		System.out.print(means.get(i)  + " ");
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
