'''
A script for plotting runtime of Sapling with different parameters across many genome lengths
Requires hard-coding the values since each one comes from its own experiment
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import math

sns.set()

# Genome lengths
xs = [4641652, 100286401, 230481012, 782475302, 2934876451, 14271578887]
logxs = [math.log(i, 2) for i in xs]

# Runtimes
sapling00001 = [45.2589, 100.803, 100.316, 137.476, 168.62, 222.989]
saplingOnePercent = [22.4503, 70.4513, 71.3830, 82.0073, 114.502, 179.87]
#saplingTenPercent = [16.4947, 62.9823, 55.3909, 94.6657, 135.960, 123.117]
saplingTwentyfivePercent = [15.7859, 45.6998, 69.0131, 62.7804, 87.2305, 115.208]
binSearch = [97.197, 171.466, 181.537, 235.254, 301.056, 369.465]
runtimeNames = ['Binary Search', 'Sapling 0.01% Overhead', 'Sapling 1% Overhead', 'Sapling 25% Overhead']

sns.lineplot(logxs, binSearch, marker = "o")
sns.lineplot(logxs, sapling00001, marker = "o")
sns.lineplot(x = logxs, y = saplingOnePercent, marker = "o")
sns.lineplot(logxs, saplingTwentyfivePercent, marker = "o")

plt.xlabel("$log_{2}(genome size)$")
plt.ylabel("$Runtime (seconds)$")
plt.title('Runtime Scaling with Genome Size')
plt.legend(labels = runtimeNames)
plt.savefig('genomelengthruntime.png')
