import pyximport; pyximport.install()
import dtw as DTW
import numpy as np


# it casts things internally, so it should be okay - should I add warnings if things get truncated?

# mess up C, send as list, then try where its not cast to a float
print("Test Bad Cost Matrix (type in array / not an np array) - expect it to handle this fine")
DTW.DTW_Cost_To_AccumCostAndSteps([[1,2,3],[4,5,6]], {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64) ,'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]]), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64) ,'SubSequence':True})

# mess up dn
print("Test bad dn (type in array / not an np array - expect it to handle this fine")
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':[1,0,1], 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64) ,'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1]), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64) ,'SubSequence':True})

# mess up dm
print("Test bad dm (type in array / not an np array) - expect it to handle this fine")
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':[1,1,0], 'dw':np.array([1,1,1], dtype=np.float64) ,'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0]), 'dw':np.array([1,1,1], dtype=np.float64) ,'SubSequence':True})

# mess up dw
print("Test bad dw (type in array / not an np array) - expect it to handle this fine")
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':[1,1,1] ,'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1]) ,'SubSequence':True})

# now try messing things up that should lead to error messagees
print("Testing missing params")
print("Nothing - expect fill in default")
DTW.DTW_Cost_To_AccumCostAndSteps([[3],[5]], {}) # send no params, should fill with default and still run
# send no dn, then no dm, then no dw
print("No dn - expect fill in default")
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64) ,'SubSequence':True})
print("No dm - expect fill in default")
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64) ,'SubSequence':True})
print("No dw - expect fill in default")
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'SubSequence':True})
print("No subsequence - expect fill in default")
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64)})


# different length dn, dm, or dw
print("Trying 3 tests where the length of dn, dm, or dw does not match, expect 3 failures")
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1,2], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64), 'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0,3], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64), 'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1], dtype=np.float64), 'SubSequence':True})

# dimensions of dn, dm, dw must be 1, C must be 2
print("Testing C wrong dimension, then dn, dm, and dw, expect 4 failures")
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([1,2,3,4,5,6], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64), 'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([[1,0],[1,0]], dtype=np.uint32), 'dm':np.array([1,1,0,0], dtype=np.uint32), 'dw':np.array([1,1,1,1], dtype=np.float64), 'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1,0], dtype=np.uint32), 'dm':np.array([[1,1],[0,0]], dtype=np.uint32), 'dw':np.array([1,1,1,1], dtype=np.float64), 'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1,1], dtype=np.uint32), 'dm':np.array([1,1,0,0], dtype=np.uint32), 'dw':np.array([[1,1],[1,1]], dtype=np.float64), 'SubSequence':True})

# testing extra dimensions but with a size of 1 (say 4x1, so the function should squeeze them into just a 1-d 4 long vector and it should work out without an error)
print("Testing trivial extra dimension (dim. of 1) - 1 failure (cost doesn't squeeze), rest run fine")
# IS THERE A GOOD WAY TO SQUEEZE COST JUST DOWN TO 2 DIMENSIONS? So to handle the extra singleton dimension
# case by squeezing, but to also still allow 1xn and nx1 cost matrices (although those have a trivial path?)

DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[[1,2,3],[4,5,6]]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64), 'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([[1],[0],[1]], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64), 'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([[[[1,1,0]]]], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64), 'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([[[1],[1],[1]]], dtype=np.float64), 'SubSequence':True})

# testing a 0, 0 step, it should alert us and fail
print("Testing 0,0 step, expect 1 failure")
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,0], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64), 'SubSequence':True})

# test for NaNs
print("Testing for NaNs, expect 1 warning and 3 failures with NaN message")
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[float('nan'),5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64), 'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':[float('nan'),0,1], 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64), 'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':[float('nan'),1,0], 'dw':np.array([1,1,1], dtype=np.float64), 'SubSequence':True})
DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':[float('nan'),1,1], 'SubSequence':True})



# step size being bigger than cost matrix should be covered by our padding?

# full one
[accumCost, steps] = DTW.DTW_Cost_To_AccumCostAndSteps(np.array([[1,2,3],[4,5,6]], dtype=np.float), {'dn':np.array([1,0,1], dtype=np.uint32), 'dm':np.array([1,1,0], dtype=np.uint32), 'dw':np.array([1,1,1], dtype=np.float64), 'SubSequence':True})


# now test get path
