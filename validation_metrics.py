from collections import defaultdict
from sklearn.metrics import mean_squared_error
import numpy as np
def rmNearZeroVarDescs(descs,ratio):
    colNb = descs.shape[1]
    to_keep=[]
    for j in range(0,colNb):
        descs_now = descs[:,j]
        d = defaultdict(int)
        for i in descs_now:
            d[i] += 1
        freqs=sorted(d.iteritems(), key=lambda x: x[1], reverse=True)
        if len(freqs) > 1: # there is more than one value..
            most_frequent1 = freqs[0][1]
            most_frequent2 = freqs[1][1]
            if (np.array(most_frequent1)/np.array(most_frequent2)) < ratio:
              to_keep.append(j)
    return np.array(to_keep)


def slope(pred,true):
    return ( np.sum(true*pred) / np.sum(pred*pred) )

def Rsquared0(pred,true):
      true_mean = true.mean()
      yr0 = pred * slope(pred,true)
      first_term = (true - yr0)**2 #(true - yr0)
      second_term = (true-true_mean)**2 #(true-true_mean)
      return (1 - ( np.sum(first_term) / np.sum(second_term) ) )

from scipy.stats import pearsonr 
def Rsquared(pred,true):
    true_mean = true.mean()
    pred_mean = pred.mean()
    first_term = np.sum( (true-true_mean) * (pred - pred_mean) )
    second_term = np.sqrt( np.sum( (true-true_mean)**2) * np.sum( (true-pred_mean)**2) )
    div = first_term / second_term
    div=pearsonr(pred,true)[0]
    return div

def Qsquared(pred,true):
    true_mean = true.mean()
    first_term = (np.abs(pred-true))**2
    second_term = (np.abs(true-true_mean))**2
    ret = 1 - (np.sum(first_term) / np.sum(second_term))
    return ret

def RMSE(pred,true):
    rmse = np.sqrt(mean_squared_error(true,pred))
    return rmse

class ValidationMetrics:
    def __init__(self, pred,true):
        self.RMSE = RMSE(pred,true)
        self.Qsquared = Qsquared(pred,true)
        self.Rsquared = Rsquared(pred,true)
        self.Rsquared0 = Rsquared0(pred,true)
        self.slope = slope(pred,true)

