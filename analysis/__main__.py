import os
from time import time
import numpy as np

import analysis.loading as load
import analysis.origin_search as origin


if __name__ == "__main__":
    
    time_ = []
    dir_path = '\\'.join(os.path.dirname(os.path.realpath(__file__)).split('\\')[:-1])

    print("Loading transfer summary")
    for list_transfer in load.serialize_files(dir_path):
        st = time()
        list_transfer = origin.merge_seq(list_transfer)

        print("  Starting hits' Blasting...")
        hitfile_path = os.path.join(dir_path, "hits")

        with open(os.path.join(hitfile_path, f"{list_transfer[0].ID}.txt"), 'w') as file :
            for transfer_elt in list_transfer:
                result = origin.blast_seq(transfer_elt)
                file.write(result)

        time_.append(time()-st)
    print("Average runtime per genome", round(np.average(time_)))
    print("Total Runtime", round(sum(time_)))
    