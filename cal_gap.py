import numpy as np
import pandas as pd
np.seterr(divide='ignore',invalid='ignore')

lamda = 0.5
alpha = 0.95

def cal_obj(WT):
    S = len(WT)
    N = int(S*(1-alpha)+1e-7)+1
    WT_sort = np.sort(WT)
    CVaR = np.mean(WT_sort[:N])
    obj = (1-lamda)*np.mean(WT)+lamda*CVaR
    return obj

array = np.zeros((10,2))
for T in range(3,9):
    S = 2**(T-1)
    data1 = np.float32(np.array(pd.read_excel('output'+str(T)+'.xlsx'))[:,1:])
    WT_plan = data1[0]

    WT_imp = np.array([data1[2**(T-2)-1+i,2*i:2*i+2] for i in range(2**(T-2))]).flatten()
    
    array[T,0] = cal_obj(WT_plan)
    array[T,1] = cal_obj(WT_imp)

print('结果为：')
print(1-array[:,1]/array[:,0])

