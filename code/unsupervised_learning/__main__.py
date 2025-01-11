import dataclasses
import compute_signatures
from scipy.base import BaseEstimator
from analysis.loading import HorizontalTransfer
import numpy as np
from datatypes import List

@dataclasses
class Predictor(BaseEstimator):
    metric:compute_signatures.metrics.Metric
    threshold:float
    window_size:int
    def predict(self,kmers_list:List[int])->np.ndarray:
        results = self.metric.slide_window(kmers_list,self.window_size)
        return results > self.threshold

def make_truth(HGTs:List[HorizontalTransfer],length:int,k:int) -> np.ndarray:
    """create a boolean array that is true when the associated kmer is in a HGT
    inputs:
    - HGTs:List[HorizontalTransfer]: a list of horisontal genes transfers into a single genome
    - length:int : the length of the genome
    - k:int: the length of the kmers
    return:
    - np.ndarray: a boolean vector
    """
    res = np.full((length,), False)
    for hgt in HGTs:
        res[hgt.start_position:(hgt.end_position+k)]=False
    return res
