# MIR_SU17_Service
The file 'Cython DTW' contains the necessary files to run, test, and benchmark the runtime efficiency of my DTW implementation. 

To just run the system you need only: dtw.pyx

An example run of the system would consist of the following python script:


'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
import numpy as np
import pyximport; pyximport.install()
import dtw


# setup your cost function
cost = np.array([[5,9,9,9], [9,8,2,9], [9,9,5,3]], dtype=np.float64)

# set parameters
dn = [1,1,2,1] # allowed steps along the rows
dm = [1,2,1,3] # allowed steps along the cols
dw = [1.0, 1.0, 2.0, 3.0] # weight of each step
subsequence = True # do subsequence matching
# create a dictionary that holds your parameters - you'll send this to the DTW function
parameter = {'dn': dn, 'dm': dm, 'dw': dw, 'SubSequence': subsequence}

# run the DTW algorithm
[accumCost, steps] = dtw.DTW_Cost_To_AccumCostAndSteps(cost, parameter)
[path, endCol, endCost] = dtw.DTW_GetPath(accumCost, steps, parameter)

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


To test the system against the Audio Labs implementation download the entire 'Cython DTW' and 'DTW_AudioLabs' files. Then run testCythonDTW_python.py in python followed by testCythonDTW.m in matlab.

To benchmark the runtime efficiency of the code, use PythonDTWTiming.py for my python implementation, and matlabDTWTiming for the matlab/c++ implementation. The matlab/c++ implementation is a bit faster.