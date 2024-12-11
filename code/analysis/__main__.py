import os
from time import time
import numpy as np

import analysis.loading as load
import analysis.origin_search as origin


if __name__ == "__main__":
    
    json_file = "transfer_summary_data_test.json"
    db_dir_name = "data_test"
    time_ = []

    dir_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

    json_path = os.path.join(dir_path, "output", json_file)
    db_path = os.path.join(dir_path, "input", "sequence_db", db_dir_name)

    print("Loading transfer summary")
    for transfer_summary in load.serialize_files(db_path, json_path=json_path):
        st = time()

        top_screen = {transfer_summary.strain : origin.screen_origins(db_path, transfer_summary)}

        best_hit_search = origin.sliding_window_search(top_screen, db_path, transfer_summary)

        print(best_hit_search)

        hitfile_path = os.path.join(dir_path, "hits")

        #
        #with open(os.path.join(hitfile_path, f"{list_transfer[0].ID}.txt"), 'w') as file :
        #    for transfer_elt in list_transfer:
        #        result = origin.blast_seq(transfer_elt)
        #        file.write(result)

        time_.append(time()-st)
        break
    print("Average runtime per genome", round(np.average(time_)))
    print("Total Runtime", round(sum(time_)))
    