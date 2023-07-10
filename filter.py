import multiprocessing
import sys
import json
import melting
import numpy as np
import pandas as pd
import subprocess

def rc(seq: str) -> str:
    '''
    Finds the reverse complement of the passed string

    Args:
        seq: sequence to process

    Returns:
        rev: the reverse complement of the passed sequence
    '''
    rev = ""
    for char in seq:
        if char == "G":
            rev = rev + "C"
        elif char == "C":
            rev = rev + "G"
        elif char == "A":
            rev = rev + "T"
        elif char == "T":
            rev = rev + "A"
        else:
            print("Cannot compute reverse complement. " +
                  seq + " is not a DNA sequence.")
    return rev

def filter_primers_into_dict(task):
    '''
    Iterates through each possible kmer in the foreground genome. Filters by melting temperature, foreground count, and background
    count. Filters can all be customized in the JSON file.

    Args:
        data: A JSON object containing hyperparameters
        fg_total_length: total length of the foreground genomes found in the JSON file
        bg_total_length: total length of the background genomes found in the JSON file

    Returns:
        primer_dict: dictionary mapping the primers that pass the filters to a list containing amount of foreground and background hits
    '''
    prefix, k, data, bg_total_length = task
    primer_dict = {}
    max_freq = 1/data['fragment_length']
    # for each prefix, iterate through each kmer length to search in the jellyfish files
    with open(prefix+'_'+str(k)+'mer_all.txt', 'r') as f_in:
        for line in f_in:
            # jellyfish dump outputs in the format "[kmer] [count]" so split on space to separate
            tupled = line.strip().split(" ")
            # isolate kmer
            kmer = tupled[0]
            # isolate count
            fg_count = int(tupled[1])
            if fg_count > data["min_fg_count"]:
                # calculate melting temp of the kmer
                tm = melting.temp(kmer)
                # check melting temp is within range
                if tm < data['max_tm'] and tm > data['min_tm']:
                    bg_count = 0
                    # iterate through each background genome given, querying the jellyfish files for counts and summing them
                    for bg_prefix in data['bg_prefixes']:
                        count_tuple = subprocess.check_output(['jellyfish', 'query', bg_prefix+'_'+str(k)+'mer_all.jf', kmer]).decode().strip().split(' ')
                        bg_count += int(count_tuple[1])
                    # final check if background frequency is below threshold, if yes then add to the primer dict
                    if bg_count/bg_total_length < max_freq:
                        primer_dict[kmer] = [fg_count, bg_count]
    return (prefix, primer_dict)

def combine_dicts(dict_list:list):
    primer_dict = {}
    for d in dict_list:
        primer_dict.update(d)
    return primer_dict

def make_tasks(mini, maxi, genomes, primer_dict, seq_lengths, prefixes):
    '''
    Initalizes the tasks array so the filter process can be multithreaded

    Args:
        data: A JSON object containing hyperparameters
        primer_dict: Dictionary mapping primers to the amount of their foreground and background hits

    Returns:
        tasks: an array of immutables of primers of each length in the range, genome paths, and sequence lengths for each
            foreground genome.
    '''
    print("Calculating position indices...")
    tasks = []
    # dicts = []
    for i, fg_genome in enumerate(genomes):
        for k in range(int(mini), int(maxi) + 1):
            primers_per_k = set(primer for primer in primer_dict[prefixes[i]] if len(primer) == k)
            reverses_per_k = set(rc(primer) for primer in primers_per_k)
            # dicts.append(get_positions((primers_per_k, k, data['fg_genomes'][i], data['fg_seq_lengths'][i])))
            if len(primers_per_k) > 0:
                tasks.append((primers_per_k, reverses_per_k, k, fg_genome, seq_lengths[i], prefixes[i]))
    return tasks

def get_positions(task):
    '''
    Reads through the foreground genome and records the indices at which each primer in the dictionary occurs then finds the lengths
    of the gaps in between the positions so the Gini index can be calculated.

    Args:
        task:
            primer_dict: dictionary mapping primers of a certain length k to empty lists to be filled with positions
            k: length of the primers in priemr_dict
            fg_genome: file path of the genome to read
            seq_length: length of the sequence in fg_genome

    Returns:
        get_gap_lengths(primer_dict, seq_length): dictionary mapping primers to their position indices and Gini index
    '''
    # breaks up the task
    primer_set, rev_set, k, fg_genome, seq_length, prefix = task

    primer_dict = {}
    rev_dict = {}
    lens = {}

    # opens the passed genome file to read
    with open(fg_genome, 'r') as f:
        index = 0
        prev = ''
        chrom_num = 0
        chr = "chr0"
        # iterates through each line in the file, concatenating the line before with the current line while removing any whitespace
        # so that kmers that cross lines can be counted. Searches the combined lines base by base checking if the substrings are in
        # the dictionary of possible primers, records any positions found.
        for line in f:
            # if line is gene header, skip
            if line[0] == ">":
                prev = ''
                if chrom_num != 0:
                    lens[chr] = index + k
                chrom_num += 1
                chr = "chr" + str(chrom_num)
                primer_dict[chr] = {}
                rev_dict[chr] = {}
                index = 0
                continue
            # combines stripped lines
            stripped = line.strip()
            combined = ''.join([prev[len(prev)-k+1:], stripped]).upper()
            # iterates through each base checking for primers
            for i in range(len(combined)-k+1):
                if combined[i:k+i] in primer_set:
                    if combined[i:k+i] not in primer_dict[chr]:
                        primer_dict[chr][combined[i:k+i]] = []
                    primer_dict[chr][combined[i:k+i]].append(index)
                elif combined[i:k+i] in rev_set:
                    if combined[i:k+i] not in rev_dict[chr]:
                        rev_dict[chr][combined[i:k+i]] = []
                    rev_dict[chr][combined[i:k+i]].append(index)
                index += 1
            prev = stripped
    lens[chr] = index + k
    return (prefix, primer_dict, rev_dict, lens)

def make_df(primers_revs_tuples:list, primer_dict:dict, prefixes:list):
    '''
    Takes in the haphazard dictionaries containing the primers of different lengths and combines them into a pandas DataFrame and a 
    dictionary mapping primers to a list their position indices

    Args:
        primers_with_ginis: list of dictionaries that map primers to a list of their position indices and their Gini index
        primer_dict: dictionary mapping primers to their foreground and background counts
        data: A JSON object containing hyperparameters

    Returns:
        sorted_df: A pandas DataFrame where each row contains a primer, its foreground counts, background counts, gini index, and 
            background/foreground count ratio, sorted by the bg/fg ratio from smallest to largest so the more specific primers are 
            at the top of the DataFrame
        primer_with_positions: Dictionary mapping each primer to a list of its position indices.
    '''
    primers_with_positions = {}
    primers_to_fg = {}
    primers_to_bg = {}
    revs_with_positions = {}
    chr_lens = {}
    for prefix in prefixes:
        primers_with_positions[prefix] = {}
        revs_with_positions[prefix] = {}
        chr_lens[prefix] = {}
        for primer in primer_dict[prefix]:
            if primer not in primers_to_fg:
                primers_to_fg[primer] = 0
                primers_to_bg[primer] = 0
            primers_to_fg[primer] += int(primer_dict[prefix][primer][0])
            primers_to_bg[primer] += int(primer_dict[prefix][primer][1])
    for trip in primers_revs_tuples:
        for chr in trip[1]:
            chr_lens[prefix][chr] = trip[3][chr]
            if chr not in primers_with_positions[prefix]:
                primers_with_positions[prefix][chr] = {}
            for primer in trip[1][chr]:
                primers_with_positions[trip[0]][chr][primer] = trip[1][chr][primer]
        for chr in trip[2]:
            if chr not in revs_with_positions[prefix]:
                revs_with_positions[prefix][chr] = {}
            for rev in trip[2][chr]:
                revs_with_positions[trip[0]][chr][rev] = trip[2][chr][rev]
    primers_as_list = []
    for primer in primers_to_bg:
        primers_as_list.append([primer, int(primers_to_fg[primer]), int(primers_to_bg[primer])])
    df = pd.DataFrame(primers_as_list, columns=['primer', 'fg_count', 'bg_count'])
    df['ratio'] = df.apply(lambda x: x['bg_count']/x['fg_count'], axis=1)
    sorted_df = df.sort_values(by=["ratio", "fg_count"], ascending=[True, False])

    return sorted_df, primers_with_positions, revs_with_positions, chr_lens

def main(data):
    bg_total_length = sum(data['bg_seq_lengths'])

    print("Filtering primers...")
    tasks = []
    for prefix in data["fg_prefixes"]:
        for k in range(int(data["min_primer_length"]), int(data["max_primer_length"]) + 1):
            tasks.append((prefix, k, data, bg_total_length))
    pool = multiprocessing.Pool(processes=int(data['cpus']))
    primer_dicts_list = pool.map(filter_primers_into_dict, tasks)

    primer_dict = {}
    num_primes = 0
    for tup in primer_dicts_list:
        num_primes += len(tup[1])
        if tup[0] not in primer_dict:
            primer_dict[tup[0]] = {}
        primer_dict[tup[0]].update(tup[1])
    print(str(num_primes) + " primers left after filtering")

    tasks = make_tasks(data["min_primer_length"], data["max_primer_length"], data["fg_genomes"], primer_dict, data["fg_seq_lengths"], data["fg_prefixes"])
    pool = multiprocessing.Pool(processes=int(data['cpus']))
    primers_revs_tuples = pool.map(get_positions, tasks)

    df, primers_with_positions, revs_with_positions, chr_lens = make_df(primers_revs_tuples, primer_dict, data['fg_prefixes'])

    if data["write"]:
        df.to_csv(data['data_dir'] + "/primers_df.csv")
        # for i, prefix in enumerate(data['fg_prefixes']):
        #     name = prefix.split('/')[-1]
        #     if len(primers_with_positions[prefix]) > 0:
        #         with open(data['data_dir'] + "/" + name + "_primers_with_positions.csv", 'w+') as f:
        #             for chr in primers_with_positions[prefix]:
        #                 f.write(">" + str(chr) + "\n")
        #                 for primer in primers_with_positions[prefix][chr]:
        #                     to_write = str(primer) + ":" + str(primers_with_positions[prefix][chr][primer]) + "\n"
        #                     f.write(to_write)
        #         with open(data['data_dir'] + "/" + name + "_reverses_with_positions.csv", 'w+') as f:
        #             for chr in revs_with_positions[prefix]:
        #                 f.write(">" + str(chr) + "\n")
        #                 for primer in revs_with_positions[prefix][chr]:
        #                     to_write = str(primer) + ":" + str(revs_with_positions[prefix][chr][primer]) + "\n"
        #                     f.write(to_write)

    return df, primers_with_positions, revs_with_positions, chr_lens

if __name__ == "__main__":
    # in_json = sys.argv[1]
    in_json = '/Users/Kaleb/Desktop/Bailey_Lab/code/newswga/params/test_params.json'
    with open(in_json, 'r') as f:
        data = json.load(f)
    main(data)
        