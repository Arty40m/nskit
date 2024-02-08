from collections import namedtuple
import math
import numpy as np

from ..parse_na import NA
from ..containers.nucleic_acid import NucleicAcid



ROUND_VALUE = 6


def _prepare_complementary_pairs(true, pred):
    if len(true)!=len(pred):
        raise ValueError((f"Compared nucleic acids must be the same length, "
                          f"got {len(true)} and {len(pred)}."))
        
    true = NA(true)
    if len(true.pairs)==0:
        raise ValueError(f"True structure has no complementary bonds.")
    
    pred = NA(pred)
    
    return true, pred


ConfusionMatrix = namedtuple("ConfusionMatrix", ["cm", "P", "N"])

def confusion_matrix(true, pred):
    true, pred = _prepare_complementary_pairs(true, pred)
    positive = set(true.pairs) # true positive bonds
    pred = set(pred.pairs)
    
    N = len(true)
    N = int((N**2 - N)/2) - len(positive) # all possible bonds minus existing
    
    TP = positive & pred # true positive (correctly predicted bonds)
    FP = pred - TP       # false positive (incorrectly predicted bonds)
    FN = positive - pred # false negative (not predicted bonds)
    TN = N - len(FP)
    
    cm = np.array([[TN, len(FP)], 
                   [len(FN), len(TP)]], 
                  dtype=np.int32)
    
    return ConfusionMatrix(cm=cm, P=positive, N=N)
    
    
###  Metrics

metrics_list = ["recall", "precision", "f1", "accuracy", "specificity", "tpr", "tnr", "fpr", "fnr"]
BinaryMetrics = namedtuple("BinaryMetrics", metrics_list)

def binary_eval(true, pred):
    conf_matrix = confusion_matrix(true, pred)
    
    tpr = _tpr(conf_matrix)
    tnr = _tnr(conf_matrix)
    return BinaryMetrics(recall=tpr, 
                         precision=_precision(conf_matrix), 
                         f1=_f_score(conf_matrix, beta=1), 
                         accuracy=_accuracy(conf_matrix), 
                         specificity=tnr, 
                         tpr=tpr, 
                         tnr=tnr, 
                         fpr=_fpr(conf_matrix), 
                         fnr=_fnr(conf_matrix)
                        )
    

### TRUE RATES

# recall - true positive rate
def _tpr(conf_matrix):
    tn, fp, fn, tp = conf_matrix.cm.ravel()
    x = tp/len(conf_matrix.P)
    return round(x, ROUND_VALUE)
    

def tpr(true, pred):
    conf_matrix = confusion_matrix(true, pred)
    return _tpr(conf_matrix)
    
recall = tpr

# specificity - true negative rate
def _tnr(conf_matrix):
    x = 1 - _fpr(conf_matrix)
    return round(x, ROUND_VALUE)
    

def tnr(true, pred):
    conf_matrix = confusion_matrix(true, pred)
    return _tnr(conf_matrix)
    
specificity = tnr
    
### FALSE RATES

# false positive rate
def _fpr(conf_matrix):
    tn, fp, fn, tp = conf_matrix.cm.ravel()
    x = fp/(fp + tn)
    return round(x, ROUND_VALUE)
    

def fpr(true, pred):
    conf_matrix = confusion_matrix(true, pred)
    return _fpr(conf_matrix)

# false negative rate
def _fnr(conf_matrix):
    x = 1 - _tpr(conf_matrix)
    return round(x, ROUND_VALUE)
    

def fnr(true, pred):
    conf_matrix = confusion_matrix(true, pred)
    return _fnr(conf_matrix)
    
    
### OTHER

# precision
def _precision(conf_matrix):
    tn, fp, fn, tp = conf_matrix.cm.ravel()
    denom = (tp + fp)
    if denom==0:
        return 0
    x = tp/denom
    return round(x, ROUND_VALUE)
    
    
def precision(true, pred):
    conf_matrix = confusion_matrix(true, pred)
    return _precision(conf_matrix)

    
# F-score
def _f_score(conf_matrix, beta):
    tn, fp, fn, tp = conf_matrix.cm.ravel()
    x = (1 + beta**2)*tp / ((1 + beta**2)*tp + (beta**2)*fn + fp)
    return round(x, ROUND_VALUE)
    

def f_score(true, pred, beta: float = 1):
    conf_matrix = confusion_matrix(true, pred)
    return _f_score(conf_matrix, beta)
    
    
# accuracy
def _accuracy(conf_matrix):
    tn, fp, fn, tp = conf_matrix.cm.ravel()
    x = (tp + tn)/(tp + tn + fp + fn)
    return round(x, ROUND_VALUE)
    

def accuracy(true, pred):
    conf_matrix = confusion_matrix(true, pred)
    return _accuracy(conf_matrix)


# Matthews correlation coefficient
def _mcc(conf_matrix):
    tn, fp, fn, tp = conf_matrix.cm.ravel()
    denom = math.sqrt( (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) )
    if denom==0:
        return 0
    x = (tp*tn + fp*fn) / denom
    return round(x, ROUND_VALUE)
    
    
def MCC(true, pred):
    conf_matrix = confusion_matrix(true, pred)
    return _mcc(conf_matrix)
    
    
    
    
    
    
    
    
    
    
    
    
    