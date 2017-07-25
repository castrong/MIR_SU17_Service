import pyximport; pyximport.install()
import DTW_Cost_To_AccumCostAndSteps as DTW
import numpy as np

DTW.DTW_Cost_To_AccumCostAndSteps([[1,2,3],[4,5,6]], {'dn':[1,0,1], 'dm':[1,1,0], 'dw':[1,1,1] ,'SubSequence':True})