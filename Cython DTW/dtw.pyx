import numpy as np
cimport numpy as np
import sys


def DTW_Cost_To_AccumCostAndSteps(C, parameter):
    '''
    Inputs
        C: The cost Matrix
    '''
    C = np.array(C)
    # make sure dn, dm, and dw are setup
    if ('dn'  in parameter.keys()):
        dn = parameter['dn']
    else:
        dn = [1, 1, 0]
    
    if 'dm'  in parameter.keys():
        dm = parameter['dm']
    else:
        dm = [1, 0, 1]
    
    if 'dw'  in parameter.keys():
        dw = parameter['dw']
    else:
        dw = [1, 1, 1]

    # add better guards, make sure C is okay / check to make sure dn / dm / dw are okay
    
    #print('dn is %s\ndm is %s\ndw is %s\n'%(str(dn), str(dm), str(dw)))

    # create matrices to store our results (D and E)
    numRows = C.shape[0] # only works with np arrays, use np.shape(x) will work on lists? want to force to use np though?
    numCols = C.shape[1]
    numDifSteps = np.size(dw)

    maxRowStep = max(dn)
    maxColStep = max(dm)

    steps = np.zeros((numRows,numCols), dtype=int)
    accumCost = np.ones((maxRowStep + numRows, maxColStep + numCols)) * float('inf')

    # essentially allow us to hop on the bottom anywhere (so could start partway through one of the signals)
    if parameter['SubSequence']:
        for col in range(numCols):
            accumCost[maxRowStep, col + maxColStep] = C[0, col]
    else:
        accumCost[maxRowStep, maxColStep] = C[0,0]

    np.set_printoptions(threshold='nan')

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
        curCol = np.argmin(accumCost[-1, :])
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





