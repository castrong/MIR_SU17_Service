import numpy as np
cimport numpy as np
cimport cython

import sys
import time

DTYPE_INT16 = np.int16
ctypedef np.int16_t DTYPE_INT16_t

DTYPE_INT32 = np.int32
ctypedef np.int32_t DTYPE_INT32_t

DTYPE_FLOAT = np.float64
ctypedef np.float64_t DTYPE_FLOAT_t

# careful, without bounds checking can mess up memory - also can't use negative indices I think (like x[-1])
@cython.boundscheck(False) # turn off bounds-checking for entire function
def DTW_Cost_To_AccumCostAndSteps(np.ndarray[DTYPE_FLOAT_t, ndim=2] C, parameter):
    '''
    Inputs
        C: The cost Matrix
    '''

    cdef np.ndarray[DTYPE_INT16_t, ndim=1] dn
    cdef np.ndarray[DTYPE_INT16_t, ndim=1] dm
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] dw
    # make sure dn, dm, and dw are setup
    if ('dn'  in parameter.keys()):
        dn = parameter['dn']
    else:
        dn = np.array([1, 1, 0], dtype=DTYPE_INT16)
    
    if 'dm'  in parameter.keys():
        dm = parameter['dm']
    else:
        dm = np.array([1, 0, 1], dtype=DTYPE_INT16)
    if 'dw'  in parameter.keys():
        dw = parameter['dw']
    else:
        dw = np.array([1, 1, 1], dtype=DTYPE_FLOAT)

    # add better guards, make sure C is okay / check to make sure dn / dm / dw are okay
    
    #print('dn is %s\ndm is %s\ndw is %s\n'%(str(dn), str(dm), str(dw)))

    # create matrices to store our results (D and E)
    cdef DTYPE_INT32_t numRows = C.shape[0] # only works with np arrays, use np.shape(x) will work on lists? want to force to use np though?
    cdef DTYPE_INT32_t numCols = C.shape[1]
    cdef DTYPE_INT16_t numDifSteps = np.size(dw)

    cdef DTYPE_INT32_t maxRowStep = max(dn)
    cdef DTYPE_INT32_t maxColStep = max(dm)

    cdef np.ndarray[DTYPE_INT16_t, ndim=2] steps = np.zeros((numRows,numCols), dtype=DTYPE_INT16)
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=2] accumCost = np.ones((maxRowStep + numRows, maxColStep + numCols), dtype=DTYPE_FLOAT) * float('inf')

    cdef DTYPE_FLOAT_t bestCost
    cdef DTYPE_INT16_t bestCostIndex
    cdef DTYPE_FLOAT_t costForStep
    cdef DTYPE_INT32_t row, col
    cdef DTYPE_INT16_t stepIndex

    # essentially allow us to hop on the bottom anywhere (so could start partway through one of the signals)
    if parameter['SubSequence']:
        for col in range(numCols):
            accumCost[maxRowStep, col + maxColStep] = C[0, col]
    else:
        accumCost[maxRowStep, maxColStep] = C[0,0]

    #np.set_printoptions(threshold='nan')

    startLoopTime = time.time()

    for row in range(maxRowStep, numRows + maxRowStep, 1):
        for col in range(maxColStep, numCols + maxColStep, 1):
            bestCost = accumCost[row, col] # initialize with what's there - so if is an entry point, then can start low
            bestCostIndex = 0
            # go through each step, find the best one
            for stepIndex in range(numDifSteps):
                costForStep = accumCost[row - dn[stepIndex], col - dm[stepIndex]] + dw[stepIndex] * C[row - maxRowStep, col - maxColStep]
                if costForStep < bestCost:
                    bestCost = costForStep
                    bestCostIndex = stepIndex
            # save the best cost and best cost index
            accumCost[row, col] = bestCost
            steps[row - maxRowStep, col - maxColStep] = bestCostIndex

    endLoopTime = time.time()

    return [accumCost[maxRowStep:, maxColStep:], steps]

def DTW_GetPath(accumCost, stepsForCost, parameter):
    '''

    Parameter should have: 'dn', 'dm', 'dw', 'SubSequence'
    '''
    numRows = accumCost.shape[0]
    numCols = accumCost.shape[1]

    # either start at the far corner (non sub-sequence)
    # or start at the lowest cost entry in the last row (sub-sequence)
    # where all of the signal along the row has been used, but only a 
    # sub-sequence of the signal along the columns has to be used
    curRow = numRows - 1
    if parameter['SubSequence']:
        curCol = np.argmin(accumCost[numCols - 1, :])
    else:
        curCol = numCols - 1

    endCol = curCol
    endCost = accumCost[curRow, curCol]

    path = np.array([[curRow],[curCol]])
    
    stepIndex = 0
    done = (parameter['SubSequence'] and curRow == 0) or (curRow == 0 and curCol == 0)
    while not done:
        if accumCost[curRow, curCol] == float('inf'):
            print('A path is not possible')

        # you're done if you've made it to the bottom left (non sub-sequence)
        # or just the bottom (sub-sequence)
        # find the step size
        curStepIndex = stepsForCost[curRow, curCol]
        curRowStep = parameter['dn'][curStepIndex]
        curColStep = parameter['dm'][curStepIndex]
        # backtrack by 1 step
        curRow = curRow - curRowStep
        curCol = curCol - curColStep
        # add your new location onto the path
        path = np.append(path, np.array([[curRow],[curCol]]), axis=1)
        # check to see if you're done
        done = (parameter['SubSequence'] and curRow == 0) or (curRow == 0 and curCol == 0)

    # reverse the path (a matrix with two rows) and return it
    return [np.fliplr(path), endCol, endCost]





