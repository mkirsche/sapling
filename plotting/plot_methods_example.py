import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
import operator

chars = ['A', 'C', 'G', 'T']

# Excerpt taken from yeast genome
s = 'AGAAGATCAAAAAAGTGAAAATTCGGTACCTTCTAAGGTTAATATGGTGAATCGCACCGATATACTGACTACGATCAAGTCATTGTCATGGCTTGACTTGATGTTGCCATTTACTATAATTCTCTCCATAATCATTGCAGTAATAATTTCTGTCTATGTGCCTTCTTCCCGTCA'

n = len(s)

ss = [s[i:] for i in range(0, n)]
ss = sorted(ss)
sa = [0 for i in range(0, n)]
for i in range(0, n):
    sa[n - len(ss[i])] = i

k = 10

bins = 8

max = 1 << (2*k)

points = []

for i in range(0, n - k + 1):
    kmer = 0
    for j in range(0, k):
        kmer *= 4
        c = s[i+j]
        if c == 'A':
            kmer += 0
        if c == 'C':
            kmer += 1
        if c == 'G':
            kmer += 2
        if c == 'T':
            kmer += 3
    points.append((kmer, sa[i]))
points.sort(key=lambda x: x[1])
xs = [x[0] for x in points]
ys = [x[1] for x in points]

for i in range(1, bins):
    plt.axvline(x=float(max) / float(bins) * i, color = 'cyan', linestyle = '--')

saplingxs = []
saplingys = []
    
for i in range(0, bins):
    minx = float(max) / float(bins) * i
    for j in range(0, len(xs)):
        if xs[j] >= minx:
            saplingxs.append(xs[j])
            saplingys.append(ys[j])
            break
saplingxs.append(xs[len(xs)-1])
saplingys.append(ys[len(ys)-1])

plt.scatter(xs, ys, color = 'pink')
plt.plot(saplingxs, saplingys, linewidth = 3.0, marker='o', color='purple')
plt.xlabel(str(k) + '-mer code')
plt.ylabel('Position in suffix array')

axes = plt.gca()
axes.set_xlim([0,max])
plt.show()
plt.savefig('methods_example.png')
            
