import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import sys

time = 'Piecewise linear time: '
buckets = 'Buckets: '

sapcolor = 'purple'
bscolor = 'pink'
bowtiecolor = 'cyan'

bstime = 47
bowtietime = 75

fn = sys.argv[1]
output = sys.argv[2]
with open(fn) as f:
    bucketlist = []
    timelist = []
    lines = f.readlines()
    for line in lines:
        if line.startswith(time):
            timelist.append(float(line[len(time):]))
        if line.startswith(buckets):
            bucketlist.append(float(line[len(buckets):]))
            
plt.scatter(bucketlist, timelist, c= sapcolor)
plt.title('Runtime of SAPLING based on sampling frequency')
plt.xlabel('log(Number of segments)')
plt.ylabel('Time (seconds per 5 million queries)')

plt.axhline(y=bstime, color=bscolor, linestyle='-')
plt.axhline(y=bowtietime, color=bowtiecolor, linestyle='-')

sappatch = mpatches.Patch(color=sapcolor, label='SAPLING')
bspatch = mpatches.Patch(color=bscolor, label='Binary search')
bowtiepatch = mpatches.Patch(color=bowtiecolor, label='Bowtie')
plt.legend(handles=[sappatch, bspatch, bowtiepatch])
plt.savefig(output)
            
