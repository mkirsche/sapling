# Makes a scatterplot from a file with an x-y pair on each line
import matplotlib.pyplot as plt
import sys
fn = sys.argv[1];
with open(fn) as f:
    lines = f.readlines()
    xs = []
    ys = []
    for line in lines:
        xs.append(int(line.split()[0]))
        ys.append(int(line.split()[1]))
    plt.scatter(xs, ys)
    plt.title('K-mers in the Human Genome')
    plt.xlabel('Kmer Code')
    plt.ylabel('Position in Suffix Array')
    plt.savefig('human.png')
