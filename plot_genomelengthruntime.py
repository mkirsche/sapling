# Makes a scatterplot from a file with an x-y pair on each line
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math
import matplotlib.patches as mpatches

sapcolor = 'pink'
bscolor = 'purple'

bsfn = sys.argv[1]
saplingfn = sys.argv[2]
bsxs = []
bsys = []

length = 'Length: '
bstime = 'Binary search time: '

with open(bsfn) as f:
    lines = f.readlines()
    lastx = 0
    for line in lines:
        if line.startswith(length):
            lastx = math.log(float(line[len(length):]), 2)
        if line.startswith(bstime):
            bsxs.append(lastx)
            bsys.append(1000*float(line[len(bstime):]))
time = 'Piecewise linear time: '

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

allxs = []
allys = []
colors = []
for i in range(0, len(bsxs)):
    allxs.append(bsxs[i])
    allys.append(bsys[i])
    colors.append(bscolor)
    
for i in range(0, len(saplingxs)):
    allxs.append(saplingxs[i])
    allys.append(saplingys[i])
    colors.append(sapcolor)
    
plt.scatter(x = allxs, y = allys, c = colors)

plt.title('Genome Length Runtime')
plt.xlabel('log(Genome Length)')
plt.ylabel('Milliseconds per 5 million queries')

sappatch = mpatches.Patch(color=sapcolor, label='SAPLING')
bspatch = mpatches.Patch(color=bscolor, label='Binary search')
plt.legend(handles=[sappatch, bspatch])

plt.savefig(sys.argv[3])
