import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math

fn = sys.argv[1]

time = 'Piecewise linear time: '
length = 'Query Length: '
saplingxs = []
saplingys = []
with open(fn) as f:
    lines = f.readlines()
    lastx = 0
    for line in lines:
        if line.startswith(length):
            lastx = float(line[len(length):])
        if line.startswith(time):
            saplingxs.append(lastx)
            saplingys.append(1000*float(line[len(time):]))
            
plt.scatter(saplingxs, saplingys)
plt.title('Query Length Runtime')
plt.xlabel('Query Length')
plt.ylabel('Milliseconds per 5 million queries')
plt.savefig(sys.argv[2])
