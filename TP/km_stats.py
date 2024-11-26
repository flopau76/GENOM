import pandas as pd
import numpy as np

from typing import List, Tuple
    
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

def load_as_matrix(jac_list : List[Tuple]) -> pd.DataFrame:
    df = pd.DataFrame(columns = sorted(list(set([jac_list[0][0]] + [elt[1] for elt in jac_list]))),
                      index = sorted(list(set([jac_list[-1][1]] + [elt[0] for elt in jac_list]))))
    
    df = df.reindex(sorted(df.columns, key=lambda x : int(x.split('_')[-1])), axis=1)
    for elt in jac_list:
        df.loc[elt[0], elt[1]] = elt[2]
    
    df = populate_matrix(df).describe()
    return df
