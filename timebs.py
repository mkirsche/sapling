import matplotlib
matplotlib.use('Agg')
import numpy as np
import random
import time
import matplotlib.pyplot as plt
import sys

ns = []
logns = []
maxn = int(sys.argv[1])
outfile = 'binsearch_runtime_' + str(maxn) + '.png'

for i in range(3, maxn):
    ns.append(2**i)
    
def binarysearch(a, want):
    lo = 0
    hi = len(a)
    while(lo < hi - 1):
        mid = (lo + hi)/2
        if(a[mid] <= want):
            lo = mid
        else:
            hi = mid
times = []
iters = 1000
for n in ns:
    a = np.arange(n)
    totTime = 0
    for i in range(0, iters):
    
        want = random.randint(0, n-1)
        start = time.time()
        x = binarysearch(a, want)
        end = time.time()
        diff = end - start
        totTime += diff
    avgTime = float(totTime) / float(iters)
    times.append(avgTime)
    logns.append(np.log2(float(n)))
plt.scatter(logns, times)
plt.ylim(0, times[len(times)-1])
plt.xlabel('Array size (log(n))')
plt.ylabel('Runtime')
plt.title('Binary Search Runtime')
plt.savefig(outfile)
