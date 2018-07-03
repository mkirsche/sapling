# Makes a scatterplot from a file with an x-y pair on each line
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

fn = sys.argv[1]
with open(fn) as f:
    lines = f.readlines()
    xs = []
    ys = []
    for line in lines:
        xs.append(float(line.split()[0]))
        ys.append(float(line.split()[1]))
    plt.scatter(xs, ys)
    plt.title('Binary Search Runtime')
    plt.xlabel('log2(genome length)')
    plt.ylabel('Milliseconds per 5 million queries')
    plt.savefig(sys.argv[2])
