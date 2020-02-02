"""
File to test models. Training is done in fit.py, and the preprocessing of SA is done in preprocess_SA.py.
This code uses PyTorch 1.1.0.
"""

import torch
import torch.utils.data as utils_data
from torch.autograd import Variable
import numpy as np
from torch import nn, optim
import matplotlib.pyplot as plt
import argparse
from collections import OrderedDict
import os
import seaborn as sns
torch.multiprocessing.set_start_method("spawn")

"""
USAGE:
python3 test.py -s <layer size> -l <number of layers> -c <convergence window> -x <total number of chunks> -y < this model's chunk> -d < location/to/save/model+plots/ > -f < location/of/preprocessed/suffixarray/ > 
e.g. python3 test.py -s 8 -l 2 -c 20 -x 10 -y 1 -d /Users/user/NN/ -f /Users/user/SuffixArray/ 

Arguments:
-s: The size of each layer in the model
-l: The number of layers in the model
-c: The convergence window to be used
-t: The tolerance used (if improvement > tolerance in the convergence window is not seen, training ends)

-x: The total number of models being used on this data
-y: Which chunk of the data test model being tested is trained on

-d: The output directory where the results and plots will be saved (must end with '/')
-f: The folder containing the preprocessed suffix array (must end with '/')

-p: If the generated plots should be shown (or just saved)
-v: Verbose flag - if all output should be printed
"""

parser = argparse.ArgumentParser()
# Model Parameters
parser.add_argument("-s", required=False, default = 20,help="Layer Size")
parser.add_argument("-l", required=False, default = 1, help="Hidden Layers")
parser.add_argument("-c", required=False, default = 10, help="Convergence Window")
parser.add_argument("-e", required=False, default = 200, help="Max Epochs")
parser.add_argument("-t", required=False, default = 0.1, help="Convergence Threshold")

parser.add_argument("-x", required=False, default = 1, help ="Number of chunks to split the data into.")
parser.add_argument("-y", required=False, default = 1, help ="Which chunk this model is for (1,..., number of chunks).")

# I/O
parser.add_argument("-d", required=True, help="Output Directory")
parser.add_argument("-f", required=True, help ="Folder containing the preprocessed suffix array.")

# Output options
parser.add_argument("-p", required=False, default = 0, help ="Show plots - 1 if yes, 0 if no.")
parser.add_argument("-v", required=False, default = 0, help ="Print full output - 1 if yes, 0 if no.")

args = parser.parse_args()

### INITIALIZE ###

# Information on layer size, number of layers, convergence and number of epochs
layerSize = int(args.s)
hiddenLayers = int(args.l)
convergenceWindow = int(args.c)
num_epochs = int(args.e)
convergenceThreshold = float(args.t)

# numChunks = total number of chinks, currentChunkNum = which chunk this model is training on.
numChunks = int(args.x)
currentChunkNum = int(args.y)

# I/O
plotDir = str(args.d)
dataDir = str(args.f)

# printing options
plotFlag = int(args.p)
print(plotFlag)
verbose = int(args.v)
if verbose == 0:
	v = False
else:
	v = True

### CUDA CHECK ###
device = torch.device("cpu")
print("Using "+ str(device))

# Print summary of the net
print("Testing a net with " + str(hiddenLayers) + " hidden layers of size " + str(layerSize) + " with a window of " + str(convergenceWindow) + ", a threshold of " + str(convergenceThreshold) + " and " +str(num_epochs) + " max epochs.")
# Print the "chunk status"
print("This model is trained on chunk " + str(currentChunkNum) + " out of " + str(numChunks))

# The filename that will be used when saving this model
outputName = "s"+str(layerSize)+"_l"+str(hiddenLayers)+"_c"+str(convergenceWindow)+"_e"+str(num_epochs)+"_t"+str(convergenceThreshold)+"_x"+str(numChunks)+"_y"+str(currentChunkNum)
print("Model filename: " + str(outputName))

# Locate the model directory
specificPlotDir = plotDir + outputName + "/"
if not os.path.exists(specificPlotDir):
	os.mkdir(specificPlotDir)

# Load the model
modelDir = specificPlotDir+outputName+'.pkl'
print("Model location: " + str(modelDir))

model = torch.load(modelDir, map_location=device)

### DATA LOADING - Load the saved numpy arrays from input folder ###

# X is the scaled decimal representations of the k-mers
Xfile = dataDir + 'x.npy'
X = np.load(Xfile)

# Y is the scaled suffix array positions
Yfile = dataDir + 'y.npy'
Y = np.load(Yfile)

# Res contains the normalized residual values
resfile = dataDir + 'res.npy'
res = np.load(resfile)

# True_res contains the original residual values (in terms of suffix array rows)
true_resfile = dataDir + 'true_res.npy'
true_res = np.load(true_resfile)

# Store some of the important values from the res array, so that we can convert error back into rows
res_min = np.min(true_res)
res_ptp = np.ptp(true_res)
min_pos = np.argmin(true_res)

### CHUNKING - by default, will do nothing to the arrays ###
print("Data loaded, getting correct chunk...")

# Divide the input data into chunks, and get appropriate chunk for this model
totalSize = len(X)
chunkSize = int(totalSize/numChunks)

# Getting the appropriate chunk
if currentChunkNum != numChunks:
	start = (currentChunkNum-1)*chunkSize
	end = (currentChunkNum)*chunkSize
	X_chunk = X[start:end]
	Y_chunk = Y[start:end]
	res_chunk =res[start:end]
	true_res_chunk = true_res[start:end]
else:
	start = (currentChunkNum-1)*chunkSize
	end = totalSize
	X_chunk = X[start:end]
	Y_chunk = Y[start:end]
	res_chunk =res[start:end]
	true_res_chunk = true_res[start:end]

# Store original, unscaled version of each
X_chunk_old = X_chunk
res_chunk_old = res_chunk

# Rescale so the values in this chunk are between 0 and 1
X_chunk = (X_chunk - np.min(X_chunk))/np.ptp(X_chunk)
res_chunk = (res_chunk - np.min(res_chunk))/np.ptp(res_chunk)

# Get the arrays we will use for testing
X = X_chunk
Y = Y_chunk
res = res_chunk
true_res = true_res_chunk


### TEST - test the model on the chunk ###
print("Testing model...")

# Pass the scaled X into the model, and retrieve the prediction
model_prediction = model(Variable(torch.Tensor(X.reshape(-1,1))))
model_prediction = model_prediction.data.numpy()

# Converting the model_prediction back into rows, so we can compute the error
# We have to do this in two steps, since we have scaled the data in two st4eps

# We first undo the scaling we performed above - this is why we store the old arrays
tmp = (model_prediction * np.ptp(res_chunk_old))+np.min(res_chunk_old)

# We then undo the scaling performed during preprocessing the suffix array (see preprocess_SA.py)
prediction = (res_ptp * tmp) + res_min

# Printing the number of predictions made by the model
print("Number of predictions: " + str(len(prediction)))


### RESULTS ###

# get the difference between the true residual values and the predicted ones
diff = np.abs(true_res - prediction)

# Output these metrics - mean, min, max, median, 25th and 75th percentile
mean_diff = np.mean(diff)
min_diff = np.min(diff)
max_diff = np.max(diff)
median_diff = np.median(diff)
diff_25 = np.percentile(diff, 25)
diff_75 = np.percentile(diff, 75)
print("Mean, Min, Max, Median, 25 and 75 percentile on relevant chunk:")
print(mean_diff, min_diff, max_diff, median_diff, diff_25, diff_75)

# Save those metrics to a file
outFile = specificPlotDir+outputName+'.txt'
f = open(outFile, "w")
f.write(str(mean_diff) + '\n')
f.write(str(min_diff)+ '\n')
f.write(str(max_diff)+ '\n')
f.write(str(median_diff)+ '\n')
f.write(str(diff_25)+ '\n')
f.write(str(diff_75)+ '\n')
f.close()

# Plot the histogram showing the spread of predictions
bins = np.arange(0, int(max_diff+2), 100)
diff_hist = plt.hist(diff, bins)

plt.title("Median Distance: " + str(median_diff))
plt.xlabel('Distance (rows)')
plt.ylabel('Count')

diffPlt = specificPlotDir+outputName+'_diffHistoWhole.png'
plt.savefig(diffPlt, bbox_inches='tight')
if plotFlag == 1:
	plt.show()
plt.clf()

# Plot the model's predicted values vs the values it was trained on
# Both the predictions and actual values are in their scaled forms
# X values remain scaled
plt.plot(X, model_prediction, c= 'r', label='prediction Line')
plt.plot(X, res, c='c', label='Prediction vs Residual (scaled)')
predPlt = specificPlotDir+outputName+'_predScaled.png'
plt.savefig(predPlt, bbox_inches='tight')
if plotFlag == 1:
	plt.show()
plt.clf()

# Plot the model's predicted values vs the values it was trained on
# Both the predictions and actual values are in their original forms
# X values remain scaled
plt.plot(X, prediction, c= 'r', label='prediction Line')
plt.plot(X, true_res, c='c', label='Prediction vs Residual (scaled)')
predPlt = specificPlotDir+outputName+'_pred.png'
plt.savefig(predPlt, bbox_inches='tight')
if plotFlag == 1:
	plt.show()
plt.clf()

# seaborn plot
sns.set()
plt.plot(X, prediction, c = 'g')
plt.plot(X, true_res, c = 'b')
plt.xlabel('K-mer Value')
plt.ylabel('Suffix Array Position')
predPlt = specificPlotDir+outputName+'_predSeaborn.png'
plt.savefig(predPlt, bbox_inches='tight')
if plotFlag == 1:
	plt.show()
plt.clf()