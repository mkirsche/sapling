import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import operator

s = 'CATATAGA'

n = len(s)

ss = [s[i:] for i in range(0, n)]
ss = sorted(ss)
sa = [0 for i in range(0, n)]
for i in range(0, n):
    sa[n - len(ss[i])] = i

k = 3

points = []
labels = []
labellocs = []
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
    if len(labels) == 0 or s[i:i+k] != labels[len(labels)-1]:
        labels.append(s[i:i+k])
        labellocs.append(kmer)
    points.append((kmer, sa[i]))
points.sort(key=lambda x: x[1])
xs = [x[0] for x in points]
ys = [x[1] for x in points]
plt.figure(figsize=(10,5))
plt.plot(xs, ys, linestyle='--', marker='o', color='b')
plt.xticks(labellocs, labels, rotation='vertical')
plt.xlabel(str(k) + '-mer code')
plt.ylabel('Position in suffix array')
plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig('kmer_example.png')
