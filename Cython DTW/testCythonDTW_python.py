import scipy.io as sp
import numpy as np
import random
import pyximport; pyximport.install()
import dtw

fileOne = sp.loadmat('test_Shostakovich_JazzSuite2_6_Waltz2_Chailly_CENS_41_10.mat')

# this example will only work (miss the high random other ones) with the right step sizes
c1 = np.random.randint(low=1, high=101, size=(9, 20))
c1[0,4] = 0.1
c1[1,6] = 0.1
c1[2,8] = 0.1
c1[3,9] = 0.1
c1[4,10] = 0.1
c1[6,11] = 0.1
c1[8,12] = 0.1


# setup your data
V = np.transpose(fileOne['f_CENS'])
Q = V[:,111:180]
fileTwo = sp.loadmat('test_Shostakovich_JazzSuite2_6_Waltz2_Yablonsky_CENS_41_10.mat')
R = np.transpose(fileTwo['f_CENS'])
c2 = 1 - np.dot(np.transpose(Q), R)




dn = np.array([1,1,2,1])
dm = np.array([1,2,1,3])
dw = np.array([1.0, 1.0, 2.0, 3.0])
subsequence = True
parameter = {'dn': dn, 'dm': dm, 'dw': dw, 'SubSequence': subsequence}

# run our version of the DTW and get:
# Cumulative cost matrices
# backtrace matrices
# optimal costs
# optimal offsets
[accumCostTwo, stepsTwo] = dtw.DTW_Cost_To_AccumCostAndSteps(c2, parameter)
[path, endCol, endCost] = dtw.DTW_GetPath(accumCostTwo, stepsTwo, parameter)
# save them as mats
sp.savemat('testCythonDTW_Subseq_python.mat', {'accumCost': accumCostTwo, 'stepMatrix': stepsTwo, 'path':path, 'offset':endCol, 'endCost':endCost})



# now try not subsequence
dn = np.array([1,1,2,1])
dm = np.array([1,2,1,3])
dw = np.array([1.0, 1.0, 2.0, 3.0])
subsequence = False
parameter = {'dn': dn, 'dm': dm, 'dw': dw, 'SubSequence': subsequence}

[accumCostTwo, stepsTwo] = dtw.DTW_Cost_To_AccumCostAndSteps(c2, parameter)
[path, endCol, endCost] = dtw.DTW_GetPath(accumCostTwo, stepsTwo, parameter)

# save them as mats
sp.savemat('testCythonDTW_nonSubseq_python.mat', {'accumCost': accumCostTwo, 'stepMatrix': stepsTwo, 'path':path, 'offset':endCol, 'endCost':endCost})



