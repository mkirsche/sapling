# Makes a scatterplot from a file with an x-y pair on each line
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math

bsfn = sys.argv[1]
saplingfn = sys.argv[2]
bsxs = []
bsys = []
with open(bsfn) as f:
    lines = f.readlines()
    for line in lines:
        bsxs.append(float(line.split()[0]))
        bsys.append(float(line.split()[1]))
time = 'Piecewise linear time: '
length = 'Length: '
saplingxs = []
saplingys = []
with open(saplingfn) as f:
    lines = f.readlines()
    lastx = 0
    for line in lines:
        if line.startswith(length):
            lastx = math.log(float(line[len(length):]), 2)
        if line.startswith(time):
            saplingxs.append(lastx)
            saplingys.append(1000*float(line[len(time):]))
print(saplingys)
allxs = []
allys = []
colors = []
for i in range(0, len(bsxs)):
    allxs.append(bsxs[i])
    allys.append(bsys[i])
    colors.append('purple')
    
for i in range(0, len(saplingxs)):
    allxs.append(saplingxs[i])
    allys.append(saplingys[i])
    colors.append('red')
    
plt.scatter(x = allxs, y = allys, c = colors)

plt.title('Genome Length Runtime')
plt.xlabel('log(Genome Length)')
plt.ylabel('Milliseconds per 5 million queries')
plt.savefig(sys.argv[3])
