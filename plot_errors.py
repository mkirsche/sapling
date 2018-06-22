# Plots prediction errors for SAPLING
# If subset is False, makes a histogram of log2(prediction error)
# If subset is true, makes a scatterplot of prediction error as function of position in first bucket

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.axes as mataxes
import matplotlib.cm as cm

import sys
import pandas
import numpy as np
import math
import random

subset = True

fn = sys.argv[1]
output = sys.argv[2]
if len(sys.argv) > 3 and sys.argv[3] == 'subset':
    subset = True
with open(fn) as f:
    df = pandas.read_table(fn)
    if subset:
        numBucketsToPlot = 10
        colors = cm.rainbow(np.linspace(0, 1, numBucketsToPlot))
        numbuckets = df['Bucket'].max()
        bucketsToPlot = random.sample(xrange(numbuckets), numBucketsToPlot)
        bucketContents = df[df['Bucket'] == bucketsToPlot[0]]
        minKmer = bucketContents['Kmer'].min()
        bucketContents['KmerDiff'] = bucketContents['Kmer'] - minKmer
        ax = bucketContents.plot(x = 'KmerDiff', y = 'Error', kind = 'scatter', c = colors[0])
        for i in range(1, numBucketsToPlot):
            bucketContents = df[df['Bucket'] == bucketsToPlot[i]]
            minKmer = bucketContents['Kmer'].min()
            bucketContents['KmerDiff'] = bucketContents['Kmer'] - minKmer
            bucketContents.plot(x = 'KmerDiff', y = 'Error', kind = 'scatter', ax = ax, c = colors[i])
    else:
        nrows = df.shape[0]
        logerrors = [float(0) for i in range(0, nrows)]
        for i in range(0, nrows):
            val = df['Error'][i]
            if val == 0:
                logerrors[i] = (float(0))
            if val > 0:
                logerrors[i] = (math.log(float(val), 2))
            if val < 0:
                logerrors[i] = (-math.log(float(-val), 2))
        df['LogError'] = logerrors
        print(nrows)
        df.plot(y='LogError', kind = 'hist')

if subset:
    plt.title('Errors within a bucket')
    plt.xlabel('Kmer code')
    plt.ylabel('Error')
else:               
    plt.title('Prediction Errors of 22mers')
    plt.xlabel('log2(Error)')
    plt.ylabel('Number of kmers')
plt.savefig(output)           
#plt.show()
