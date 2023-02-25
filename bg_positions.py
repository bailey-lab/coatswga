import sys
import pandas as pd
import numpy as np
import os
import json
import multiprocessing
import h5py

def get_bg_positions(task):
    primer_dict, k, bg_prefix, bg_genome = task
    print("Getting positions...")
    with open(bg_genome, 'r') as f:
        index = 0
        prev = ''
        for line in f:
            if line[0] == ">":
                prev = ''
                continue
            stripped = line.strip()
            combined = ''.join([prev[len(prev)-k+1:], stripped]).upper()
            for i in range(len(combined)-k+1):
                if combined[i:k+i] in primer_dict:
                    primer_dict[combined[i:k+i]].append(index)
                index += 1
            prev = stripped
    print("Writing to h5...")
    if os.path.exists(bg_prefix + '_' + str(k) + 'mer_positions.h5'):
        os.remove(bg_prefix + '_' + str(k) + 'mer_positions.h5')
    file = h5py.File(bg_prefix + '_' + str(k) + 'mer_positions.h5', 'a')
    for primer in primers_per_k:
        positions = primers_per_k[primer]
        if primer not in file:
            file.create_dataset(primer, data=positions)
        else:
            del file[primer]
            file.create_dataset(primer, data=positions)
    file.close()

if __name__ == "__main__":
    in_json = sys.argv[1]
    # in_json = '/Users/Kaleb/Desktop/Bailey_Lab/code/newswga/new_src/test_params.json'
    with open(in_json, 'r') as f:
        data = json.load(f)
    df = pd.read_csv(os.path.join(data['data_dir'], "step2_df.csv"))
    tasks = []
    for i, bg_prefix in enumerate(data['bg_prefixes']):
        for k in range(int(data["min_primer_length"]), int(data["max_primer_length"]) + 1):
            primers_per_k = {primer: [] for primer in df['primer'] if len(primer) == k}
            # get_bg_positions((primers_per_k, k, bg_prefix, data['bg_genomes'][i]))
            if primers_per_k != {}:
                tasks.append((primers_per_k, k, bg_prefix, data['bg_genomes'][i]))
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    pool.map(get_bg_positions, tasks)


# import sys
# import pandas as pd
# import numpy as np
# import os
# import json
# import multiprocessing
# import h5py

# def get_bg_positions(task):
#     primer_dict, bg_prefix, bg_genome = task
#     with open(bg_genome, 'r') as f:
#         index = 0
#         prev = ''
#         for line in f:
#             if line[0] == ">":
#                 prev = ''
#                 continue
#             stripped = line.strip()
#             combined = ''.join([prev[len(prev)-k+1:], stripped]).upper()
#             for i in range(len(combined)-k+1):
#                 if combined[i:k+i] in primer_dict:
#                     primer_dict[combined[i:k+i]].append(index)
#                 index += 1
#             prev = stripped
#     write_to_h5(bg_prefix, primer_dict)

# def write_to_h5(bg_prefix:str, primer_dict):
#     for k in range(int(data["min_primer_length"]), int(data["max_primer_length"]) + 1):
#         primers_per_k = {primer: primer_dict[primer] for primer in primer_dict if len(primer) == k}
#         if os.path.exists(bg_prefix + '_' + str(k) + 'mer_positions.h5'):
#             file = open(bg_prefix + '_' + str(k) + 'mer_positions.h5')
#         else:
#             file = h5py.File(bg_prefix + '_' + str(k) + 'mer_positions.h5', 'a')
#         for primer in primers_per_k:
#             positions = primers_per_k[primer]
#             if primer not in file:
#                 file.create_dataset(primer, data=positions)
#             else:
#                 del file[primer]
#                 file.create_dataset(primer, data=positions)
#         file.close()

# if __name__ == "__main__":
#     in_json = sys.argv[1]
#     # in_json = '/Users/Kaleb/Desktop/Bailey_Lab/code/newswga/new_src/test_params.json'
#     with open(in_json, 'r') as f:
#         data = json.load(f)
#     df = pd.read_csv(os.path.join(data['data_dir'], "step2_df.csv"))
#     tasks = []
#     for i, bg_prefix in enumerate(data['bg_prefixes']):
#         primer_dict = {primer: [] for primer in df['primer']}
#         get_bg_positions((primer_dict, bg_prefix, data['bg_genomes'][i]))