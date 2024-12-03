import pandas as pd
import numpy as np

import heapq, math
from typing import List, Tuple
from collections import Counter, deque
    
from TP.kmers import stream_kmers

def populate_matrix(df : pd.DataFrame) -> pd.DataFrame:
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            if i == j:
                df.iat[i,j] = np.nan
            elif str(df.iat[i,j]) == 'nan':
                df.iat[i,j] = df.iat[j, i]

    for col in df.columns:
        df[col] = df[col].astype(float)
    return df

def vect_stats(df : pd.DataFrame, df_stats : pd.DataFrame) -> pd.DataFrame:
    df = df.apply(lambda x : np.abs(x-x.mean()), axis=1)
    return df

def load_as_matrix(jac_list : List[Tuple]) -> pd.DataFrame:
    df = pd.DataFrame(columns = sorted(list(set([jac_list[0][0]] + [elt[1] for elt in jac_list]))),
                      index = sorted(list(set([jac_list[-1][1]] + [elt[0] for elt in jac_list]))))
    
    df = df.reindex(sorted(df.columns, key=lambda x : int(x.split('_')[-1])), axis=1)
    df = df.reindex(sorted(df.columns, key=lambda x : int(x.split('_')[-1])), axis=0)
    for elt in jac_list:
        df.loc[elt[0], elt[1]] = elt[2]
    
    df = populate_matrix(df)
    return df

def whole_genome_signature(kmers_list : list) -> Counter :
    return Counter(kmers_list)


def window_slider(kmers_list : list[int], signature_km_freq : dict, width : int = 2000, heap_size : int = 5) -> float:
    heap = [(0, 0) for _ in range(heap_size)]
    kmers_list = [signature_km_freq[i] for i in kmers_list]

    all_avg = -np.log10(np.convolve(kmers_list, np.ones(width), 'valid')/width)
    for i in range(width, len(all_avg)-width+1):
        heapq.heappushpop(heap, (all_avg[i], i))

    res = {key:val for val, key in heap}
    return res, all_avg
