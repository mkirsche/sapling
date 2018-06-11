import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

time = 'Piecewise linear time: '
buckets = 'Buckets: '

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
            
plt.scatter(bucketlist, timelist)
plt.title('Runtime of SAPLING based on sampling frequency')
plt.xlabel('Log_2(Number of buckets)')
plt.ylabel('Time (seconds per 5 million queries)')
plt.savefig(output)
            
