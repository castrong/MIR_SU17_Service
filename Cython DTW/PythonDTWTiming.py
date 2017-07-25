import scipy.io as sp
import numpy as np
import random
import pyximport; pyximport.install()
import dtw
import time

fileOne = sp.loadmat('test_Shostakovich_JazzSuite2_6_Waltz2_Chailly_CENS_41_10.mat')

# setup your data
V = np.transpose(fileOne['f_CENS'])
Q = V[:,111:180]
fileTwo = sp.loadmat('test_Shostakovich_JazzSuite2_6_Waltz2_Yablonsky_CENS_41_10.mat')
R = np.transpose(fileTwo['f_CENS'])
c2 = 1 - np.dot(np.transpose(Q), R)

# setup parameters

dn = np.array([1,1,2,1])
dm = np.array([1,2,1,3])
dw = np.array([1.0, 1.0, 2.0, 3.0])
subsequence = True
parameter = {'dn': dn, 'dm': dm, 'dw': dw, 'SubSequence': subsequence}

# loop through and do it a few times
numIterations = 100
firstStepTimes = []
secondStepTimes = []
totalTimes = []

for i in range(numIterations):
	# do dtw
	startTime = time.time()
	[accumCostTwo, stepsTwo] = dtw.DTW_Cost_To_AccumCostAndSteps(c2, parameter)
	firstStepEndTime = time.time()
	[path, endCol, endCost] = dtw.DTW_GetPath(accumCostTwo, stepsTwo, parameter)
	secondStepEndTime = time.time()

	firstStepTimes.append(firstStepEndTime - startTime)
	secondStepTimes.append(secondStepEndTime - firstStepEndTime)
	totalTimes.append(secondStepEndTime - startTime)

print("After %g iterations of subsequence matching"%(numIterations))
print("First Step: %f"%(np.mean(firstStepEndTime - startTime)))
print("Second Step: %f"%(np.mean(secondStepEndTime - firstStepEndTime)))
print("Total: %f"%(np.mean(secondStepEndTime - startTime)))




# setup parameters

dn = np.array([1,1,2,1])
dm = np.array([1,2,1,3])
dw = np.array([1.0, 1.0, 2.0, 3.0])
subsequence = False
parameter = {'dn': dn, 'dm': dm, 'dw': dw, 'SubSequence': subsequence}

# loop through and do it a few times
numIterations = 100
firstStepTimes = []
secondStepTimes = []
totalTimes = []

for i in range(numIterations):
	# do dtw
	startTime = time.time()
	[accumCostTwo, stepsTwo] = dtw.DTW_Cost_To_AccumCostAndSteps(c2, parameter)
	firstStepEndTime = time.time()
	[path, endCol, endCost] = dtw.DTW_GetPath(accumCostTwo, stepsTwo, parameter)
	secondStepEndTime = time.time()

	firstStepTimes.append(firstStepEndTime - startTime)
	secondStepTimes.append(secondStepEndTime - firstStepEndTime)
	totalTimes.append(secondStepEndTime - startTime)

print("After %g iterations of non-subsequence matching"%(numIterations))
print("First Step: %f"%(np.mean(firstStepEndTime - startTime)))
print("Second Step: %f"%(np.mean(secondStepEndTime - firstStepEndTime)))
print("Total: %f"%(np.mean(secondStepEndTime - startTime)))