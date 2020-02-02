"""
File to preprocess suffix array input.
Requires the processed suffix array as input.
"""

import numpy as np
import os
import argparse

"""
Example usage:
python preprocess.py -f </path/to/sa.key.row> -d </path/to/preprocessed/SA/>

Arguments:
-f: Path to the folder containing the suffix array file (.key.row; must end in '/')
-d: Path to the folder where the numpy arrays should be saved (must end in '/')
-v: Verbose flag
"""
parser = argparse.ArgumentParser()
parser.add_argument("-f", required=True, help ="Input File")
parser.add_argument("-d", required=True, help="Output Directory")
parser.add_argument("-v", required=False, default = 0, help ="Print full output - 1 is yes, 0 if no.")
args = parser.parse_args()

outDir = str(args.d)
dataFile = str(args.f)
verbose = int(args.v)
if verbose == 0:
	v = False
else:
	v = True


### DATA LOADING AND REMOVING OUTLIERS ###

# Load the data into a numpy array
data = np.loadtxt(dataFile, usecols=(0,1),skiprows=1, delimiter=" ")
data = data.T
print("Processing data...")

# Turn it into a list of lists
X = np.array([[i]for i in data[0]])
y = np.array([[i]for i in data[1]])

# Removing the outliers - these are the rows with strings that are shorter than the other k-mers
# We can identify them as values that are smaller than the value before them (as they come from a shorter string)
if v:
	print("X before filtering:")
	print(X)
	print("y before filtering:")
	print(y)
toDelete = []
if v:
	print("Filtered Values: Index, X-value, y-value")
for i in range(1,len(X)):
	if X[i][0] < X[i-1][0]:
		if v:
			print(i, X[i][0], y[i][0])
		toDelete.append(i)
toDelete = np.asarray(toDelete)
if v:
	print("Indices to be deleted")
	print(toDelete)
	print("X values to be deleted")
	print(X[toDelete])
	print("Y values to be deleted")
	print(y[toDelete])
	print("Original array length (X,y):")
	print(len(X), len(y))
X = np.delete(X, toDelete)
y = np.delete(y, toDelete)

if v:
	print("New array length (X,y):")
	print(len(X), len(y))

	print("X after filtering:")
	print(X)
	print("y after filtering:")
	print(y)

X = np.array([[i] for i in X])
y = np.array([[i] for i in y])

if v:
	print("X after filtering and reshape:")
	print(X)
	print("y after filtering and reshape:")
	print(y)

# Iterate through and make sure no such outlier exists - values should be in increasing order
for i in range(1,len(X)):
	if X[i][0] < X[i-1][0]:
		print("Outlier found!")
		print(i)

### NORMALIZING AND SCALING ###

# Scaling X to bring it between 0 and 1.
X = X/np.max(X)

# Generating the residual - done by finding the straight line between the first and last points.

# compute gradient
dX = X[len(X)-1][0] - X[0][0]
dy = y[len(y)-1][0] - y[0][0]
m = dy/dX
# get intercept
c = y[0][0] - X[0][0]*m

m_val = m
c_val = c

# Compute residual values by taking the difference between the straight line function and the real value
res = []
for idx in range(len(X)):
	lineVal = X[idx][0] * m + c
	trueVal = y[idx][0]
	resVal = lineVal - trueVal
	res.append([resVal])

# Keep a copy of the residual, since we will be scaling these values
true_res = res

# Scale the residual value - get the minimum value and the range (min to max)
res_min = np.min(res)
res_ptp = np.ptp(res)
min_pos = np.argmin(res)

# Scale by substracting the min value (making everything positive) and dividing by the range (between 0 and 1)
res = (res - np.min(res))/np.ptp(res)


### SAVE TO FILE ###
print("Creating arrays....")

xName = outDir + 'x'
yName = outDir + 'y'
resName = outDir + 'res'
true_resName = outDir + 'true_res'

np.save(xName, X)
np.save(yName, y)
np.save(resName, res)
np.save(true_resName, true_res)