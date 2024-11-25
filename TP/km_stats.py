import pandas as pd
import numpy as np

from typing import List, Tuple

def load_as_matrix(jac_list : List[Tuple]):
    df = pd.DataFrame(columns=[elt[1] for elt in jac_list],
                      index = [elt[0] for elt in jac_list])
    return df

load_as_matrix()