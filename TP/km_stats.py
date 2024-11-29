import pandas as pd
import numpy as np

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


def window_slider(kmers_list : list[int], signature_km_freq : dict, width : int = 200) -> float:
    window = deque([signature_km_freq[i] for i in kmers_list[0:width]])
    all_avg = [np.log10(sum(window)/width)]
    for i in range(width,len(kmers_list)-width+1):
        window.popleft()
        window.append(signature_km_freq[kmers_list[i]])

        all_avg.append(np.log10(sum(window)/width))

    return all_avg
