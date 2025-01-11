import dataclasses
import compute_signatures
import numpy as np
from datatypes import List
@dataclasses
class Predictor:
    metric:compute_signatures.metrics.Metric
    threshold:float
    window_size:int
    def predict(self,kmers_list:List[int])->np.ndarray:
        results = self.metric.slide_window(kmers_list,self.window_size)
        return results > self.threshold
