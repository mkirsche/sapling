'''
Code for plotting k-mer uniqueness ratio
Assumes count.cpp has already beeen run to put data in (k, # different kmers, # unique kmers, # kmers) format
'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

files = ['ecoli.kmers', 'celegans.kmers', 'chr1.kmers', 'tomato.kmers', 'human.kmers', 'wheat.kmers']
names = ['E. coli', 'C. elegans', 'H. sapiens (chr1)', 'S. lycopersicum', 'H. sapiens', 'T. aestivum']
colors = ['green', 'blue', 'purple', 'red', 'orange', 'yellow']
xs = []
ys = []
for fn in files:
  with open(fn) as f:
    curxs = []
    curys = []
    for line in f:
      tokens = line.split(' ')
      if len(tokens) == 4:
        x = int(tokens[0])
        y = float(tokens[2]) / float(tokens[3])
        if x > 300:
          continue
        curxs.append(x)
        curys.append(y)
    xs.append(curxs)
    ys.append(curys)
sns.set()
for i in range(0, len(xs)):
  sns.lineplot(xs[i], ys[i], color = colors[i])
plt.title('K-mer Uniqueness Ratio')
plt.xlabel('K-mer Length')
plt.ylabel('Proportion of Unique K-mers')
plt.legend(labels = names)
plt.savefig('uniquekmers.png')
