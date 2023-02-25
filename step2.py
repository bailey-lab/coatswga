import multiprocessing
import sys
import json
import melting
import numpy as np
import pandas as pd
import os
import subprocess
import h5py

def filter_primers_into_dict(data, fg_total_length, bg_total_length):
    primer_dict = {}
    print("Filtering primers...")
    for prefix in data["fg_prefixes"]:
        for k in range(int(data["min_primer_length"]), int(data["max_primer_length"]) + 1):
            with open(prefix+'_'+str(k)+'mer_all.txt', 'r') as f_in:
                for line in f_in:
                    tupled = line.strip().split(" ")
                    kmer = tupled[0]
                    fg_count = tupled[1]
                    fg_freq = int(fg_count)/fg_total_length
                    if fg_freq > data["min_fg_freq"]:
                        tm = melting.temp(kmer)
                        if tm < data['max_tm'] and tm > data['min_tm']:
                            bg_count = 0
                            for bg_prefix in data['bg_prefixes']:
                                count_tuple = subprocess.check_output(['jellyfish', 'query', bg_prefix+'_'+str(k)+'mer_all.jf', kmer]).decode().strip().split(' ')
                                bg_count += int(count_tuple[1])
                            if bg_count/bg_total_length < data["max_bg_freq"]:
                                primer_dict[kmer] = [fg_count, bg_count]
    return primer_dict

def calc_gini(data, primer_dict):
    print("Calculating Gini indices...")
    tasks = []
    dicts = []
    for i, fg_prefix in enumerate(data['fg_prefixes']):
        for k in range(int(data["min_primer_length"]), int(data["max_primer_length"]) + 1):
            primers_per_k = {primer: [[]] for primer in primer_dict if len(primer) == k}
            # dicts.append(get_positions((primers_per_k, k, data, fg_prefix, data['fg_genomes'][i], data['fg_seq_lengths'][i])))
            tasks.append((primers_per_k, k, data['fg_genomes'][i], data['fg_seq_lengths'][i]))
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    return pool.map(get_positions, tasks)
    # return dicts


def get_positions(task):
    primer_dict, k, fg_genome, seq_length = task
    with open(fg_genome, 'r') as f:
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
                    primer_dict[combined[i:k+i]][0].append(index)
                index += 1
            prev = stripped
    return get_gap_lengths(primer_dict, seq_length)

def get_gap_lengths(primer_dict, seq_length):
    for primer in primer_dict:
        if type(primer_dict[primer][0]) != list or primer_dict[primer] == []:
            continue
        positions = primer_dict[primer][0]
        differences = [a_i - b_i for a_i, b_i in zip(positions[1:], positions[:-1])]
        differences.append(seq_length-positions[-1])
        differences.append(positions[0])
        primer_dict[primer].append(gini_exact(differences))
    return primer_dict

def gini_exact(array):
    """Calculate the Gini coefficient of a numpy array. Based on bottom equation from
    http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    All values are treated equally, arrays must be 1d

    Args:
        array (array): One-dimensional array or list of floats or ints.

    Returns:
        gini_index: The gini index of the array given.
    """
    if len(array) == 0:
        return 0
    array = np.array(array, dtype=float)
    if np.amin(array) < 0:
        # Values cannot be negative:
        array -= np.amin(array)
    # Values cannot be 0:
    array += 0.0000001
    # Values must be sorted:
    array = np.sort(array)
    # Index per array element:
    index = np.arange(1, array.shape[0] + 1)
    # Number of array elements:
    n = array.shape[0]
    # Gini coefficient:
    return ((np.sum((2 * index - n - 1) * array)) / (n * np.sum(array)))

def make_df(primers_with_ginis:dict):
    primers_as_list = []
    for diction in primers_with_ginis:
        for primer in diction:
            primers_as_list.append([primer, int(primer_dict[primer][0]), int(primer_dict[primer][1]), diction[primer][1]])
    df = pd.DataFrame(primers_as_list, columns=['primer', 'fg_count', 'bg_count', 'gini'])
    df['gini_bool'] = df.apply(lambda x: x['gini'] is not None and x['gini'] < data['max_gini'], axis=1)
    df['ratio'] = df.apply(lambda x: x['bg_count']/x['fg_count'], axis=1)
    sorted_df = df.sort_values(by=["ratio"])[:data['max_primer']]
    return sorted_df

def make_h5(primers_with_positions:dict, df:pd.DataFrame):
    for fg_prefix in data['fg_prefixes']:
        for primers_per_k in primers_with_positions:
            if primers_per_k == {}:
                continue
            k = len(list(primers_per_k.keys())[0])
            file = h5py.File(fg_prefix + '_' + str(k) + 'mer_positions.h5', 'a')
            for primer in primers_per_k:
                positions = primers_per_k[primer][0]
                if primer not in file:
                    file.create_dataset(primer, data=positions)
                else:
                    del file[primer]
                    file.create_dataset(primer, data=positions)
            file.close()

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
    for primer in primer_dict:
        positions = primer_dict[primer]
        if primer not in file:
            file.create_dataset(primer, data=positions)
        else:
            del file[primer]
            file.create_dataset(primer, data=positions)
    file.close()

def start_bg_positions(data, df):
    tasks = []
    for i, bg_prefix in enumerate(data['bg_prefixes']):
        for k in range(int(data["min_primer_length"]), int(data["max_primer_length"]) + 1):
            primers_per_k = {primer: [] for primer in df['primer'] if len(primer) == k}
            # get_bg_positions((primers_per_k, k, bg_prefix, data['bg_genomes'][i]))
            if primers_per_k != {}:
                tasks.append((primers_per_k, k, bg_prefix, data['bg_genomes'][i]))
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    pool.map(get_bg_positions, tasks)

if __name__ == "__main__":
    in_json = sys.argv[1]
    # in_json = '/Users/Kaleb/Desktop/Bailey_Lab/code/newswga/new_src/test_params.json'
    with open(in_json, 'r') as f:
        data = json.load(f)
    fg_total_length = sum(data['fg_seq_lengths'])
    bg_total_length = sum(data['bg_seq_lengths'])
    primer_dict = filter_primers_into_dict(data, fg_total_length, bg_total_length)
    primers_with_ginis = calc_gini(data, primer_dict)
    df = make_df(primers_with_ginis)
    make_h5(primers_with_ginis, df)
    start_bg_positions(data, df)
    df.to_csv(os.path.join(data['data_dir'], "step2_df.csv"))
