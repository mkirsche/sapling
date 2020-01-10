import java.io.File;
import java.io.FileInputStream;
import java.util.HashMap;
import java.util.Scanner;

public class AlignmentQuality {
	static int MAX_ERROR = 10;
public static void main(String[] args) throws Exception
{
	String fn = "/home/mkirsche/git/SaplingEval/ecoli_mine.sam";
	String truth = "/home/mkirsche/ecoli_sim.sam";
	if(args.length >= 2)
	{
		fn = args[0];
		truth = args[1];
	}
	
	// Scan through true alignments and store the position for each read name
	Scanner input = new Scanner(new FileInputStream(new File(truth)));
	HashMap<String, String> readToChr = new HashMap<String, String>();
	HashMap<String, Integer> readToPos = new HashMap<String, Integer>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0 || line.startsWith("@"))
		{
			continue;
		}
		String[] tokens = line.split("\t");
		String readName = tokens[0];
		String chrName = tokens[2];
		int chrPos = Integer.parseInt(tokens[3]);
		readToChr.put(readName, chrName);
		readToPos.put(readName, chrPos);
	}
	input.close();
	
	// Scan through my alignments and see how many are right
	input = new Scanner(new FileInputStream(new File(fn)));
	int good = 0, bad = 0;
	int unaligned = 0;
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0 || line.startsWith("@"))
		{
			continue;
		}
		String[] tokens = line.split("\t");
		String readName = tokens[0];
		String chrName = tokens[2];
		if(chrName.equals("*"))
		{
			unaligned++;
			continue;
		}
		int chrPos = Integer.parseInt(tokens[3]);
		String trueChrName = readToChr.get(readName);
		int truePos = readToPos.get(readName);
		if(chrName.equals(trueChrName) && Math.abs(chrPos - truePos) <= MAX_ERROR)
		{
			good++;
		}
		else
		{
			bad++;
		}
	}
	System.out.println("Correctly aligned reads: " + good);
	System.out.println("Incorrectly aligned reads: " + bad);
	System.out.println("Unaligned reads: " + unaligned);
	input.close();
}
}
