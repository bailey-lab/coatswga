import multiprocessing
import sys
import json
import melting
import numpy as np
import pandas as pd
import subprocess
from time import perf_counter as pc

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
            rev = "C" + rev
        elif char == "C":
            rev = "G" + rev
        elif char == "A":
            rev = "T" + rev
        elif char == "T":
            rev = "A" + rev
        else:
            print("Cannot compute reverse complement. " +
                  seq + " is not a DNA sequence.")
    return rev

def filter_primers_into_dict(task):
    '''
    Iterates through each possible kmer in the foreground genome. Filters by melting temperature, foreground count, and background
    count. Filters can all be customized in the JSON file.

    Args:
        task:
            prefix: Foreground sequence prefix to use to open kmer files
            k: Lenth of kmers to look at
            data: Json file containing hyperparameters
            bg_total_length: Sum of background genomes lengths, used to filter by background hits

    Returns:
        tuple:
            prefix: The prefix passed to the function, used to separate kmer counts from different genomes
            primer_dict: Dictionary mapping the primers that pass the filters to a list containing amount of foreground and background hits
    '''
    prefix, k, data = task
    primer_set = set()
    primer_dict = {}
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
                    # final check if background frequency is below threshold, add to the primer dict if yes
                    if bg_count/fg_count < data['max_ratio']:
                        primer_dict[kmer] = [fg_count, bg_count]
                        primer_set.add(kmer)
                        primer_set.add(rc(kmer))
    if data['verbose']:
        print(f"{len(primer_dict)} {k}-mers in {prefix}")
    return (prefix, primer_dict, primer_set, k)

def get_positions(task):
    '''
    Reads through the foreground genome and records the indices at which each primer and its reverse complement occurs

    Args:
        task:
            primer_set: Set containing all primers and reverse complements of length k
            k: Length of the primers given
            genome: File path of the genome to read
            prefix: Prefix of the foreground genome, used to ID the immutable 

    Returns:
        immutable:
            prefix: Prefix of the foreground genome, used to ID the immutable 
            primer_dict: Dictionary of primers of length k mapping to a list of positions where that primer is found in each chromosome for the given genome
            lens: Dicitonary mapping chromosome numbers to their lengths
    '''
    # t0 = pc()
    # breaks up the task
    primer_set, k, genome, prefix = task
    # initialize variables to hold sequence 
    seq, name_list, seq_list, = '', [], []
    for line in open(genome):
        line=line.strip()
        if '>' in line:
            name_list.append(line[1:])
            if len(seq)>0:
                seq_list.append(seq)
                seq=''
        else:
            seq=seq+line.upper()
    seq_list.append(seq)
    genome_list = [[name, seq_list[name_number]] for name_number, name in enumerate(name_list)]

    kmer_dict = {}
    lens = {}
    for chrom, seq in genome_list:
        kmer_dict[chrom] = {}
        lens[chrom] = len(seq)
        for index in range(len(seq)-k+1):
            kmer = seq[index:index+k]
            if kmer in primer_set:
                if kmer not in kmer_dict[chrom]:
                    kmer_dict[chrom][kmer] = []
                kmer_dict[chrom][kmer].append(index)
    # print(f"{len(primer_set)} {k}mers found in {pc() - t0}")
    return (prefix, kmer_dict, lens)

def make_df(immut_list:list, primer_dict:dict, prefixes:list):
    '''
    Takes in the haphazard dictionaries containing the primers of different lengths and combines them into a pandas DataFrame and a 
    dictionary mapping primers to a list their position indices

    Args:
        immut_list: List of immutables outputted by get_positions()
        primer_dict: Dictionary mapping primers to their foreground and background counts
        prefixes: List of foreground prefixes

    Returns:
        sorted_df: A pandas DataFrame where each row contains a primer, its foreground counts, background counts, and background/foreground count ratio, 
                    sorted by the bg/fg ratio from smallest to largest then by fg count so the most common of the most specific primers are at the top of the DataFrame
        primer_with_positions: Dictionary mapping each primer to a list of its position indices within each chromosome in each foreground genome.
        revs_with_positions: Same as primer_with_positions but with reverse complements
        chr_lens: Dictionary mapping each chromosome of each fg prefix to its length
    '''
    # initialize data structures to hold the combined primer dicts
    primers_with_positions = {}
    primers_to_fg = {}
    primers_to_bg = {}
    chr_lens = {}
    for prefix in prefixes:
        primers_with_positions[prefix] = {}
        chr_lens[prefix] = {}
        # adds the fg count and bg count to their respective dictionaries for each primer to sum fg and bg counts across each prefix
        for primer in primer_dict[prefix]:
            if primer not in primers_to_fg:
                primers_to_fg[primer] = 0
                primers_to_bg[primer] = 0
            primers_to_fg[primer] += int(primer_dict[prefix][primer][0])
            primers_to_bg[primer] += int(primer_dict[prefix][primer][1])
    # iterates through each immutable in the list
    for trip in immut_list:
        pref = trip[0]
        # iterates through each chromosome in the dictionary of primers
        for chr in trip[1]:
            # adds the length of the chromosome to the above dictionary
            chr_lens[pref][chr] = trip[2][chr]
            # if the chromosome has not been added already, initialize an empty dictionary to fill with primers
            if chr not in primers_with_positions[pref]:
                primers_with_positions[pref][chr] = {}
            # for each primer in the chromosome add the positions to the central dictionary
            for primer in trip[1][chr]:
                primers_with_positions[pref][chr][primer] = trip[1][chr][primer]
    # initialize list to use for dataframe
    primers_as_list = []
    # since the fg and bg dictionaries contain the same primers, only need to iterate through one
    for primer in primers_to_bg:
        # for each primer add a row to the list containing the primer, fg count, and bg count
        primers_as_list.append([primer, int(primers_to_fg[primer]), int(primers_to_bg[primer])])
    # make a DataFrame from the list
    df = pd.DataFrame(primers_as_list, columns=['primer', 'fg_count', 'bg_count'])
    # calculate the bg/fg ratio for each primer
    df['ratio'] = df.apply(lambda x: x['bg_count']/x['fg_count'], axis=1)
    # sort the df first by ratio (low to high) then by fg count (high to low) so the most common and most specific primers will be at the top
    # sorted_df = df.sort_values(by=["ratio", "fg_count"], ascending=[True, False])
    df['sort_val'] = df.apply(lambda x: x['fg_count']/x['ratio'])
    sorted_df = df.sort_values(by='sort_val', ascending=True)

    return sorted_df, primers_with_positions, chr_lens

def main(data):
    t0 = pc()

    print("Filtering primers...")
    tasks = []
    for prefix in data["fg_prefixes"]:
        for k in range(int(data["min_primer_length"]), int(data["max_primer_length"]) + 1):
            tasks.append((prefix, k, data))
    pool = multiprocessing.Pool(processes=int(data['cpus']))
    primer_dicts_list = pool.map(filter_primers_into_dict, tasks)

    primer_dict = {}
    primer_sets = {}
    num_primes = 0
    for tup in primer_dicts_list:
        num_primes += len(tup[1])
        if tup[0] not in primer_dict:
            primer_dict[tup[0]] = {}
            primer_sets[tup[0]] = {}
        primer_dict[tup[0]].update(tup[1])
        primer_sets[tup[0]][tup[3]] = tup[2]
    print(str(num_primes) + " primers left after filtering")
    
    t1 = pc()
    if data['verbose']:
        print("Time to filter primers: " + str(round(t1 - t0, 4)))

    print("Calculating position indices...")
    ta = pc()
    pool = multiprocessing.Pool(processes=int(data['cpus']))
    tasks = []
    for i, fg_genome in enumerate(data['fg_genomes']):
        prefix = data['fg_prefixes'][i]
        for k in range(int(data['min_primer_length']), int(data['max_primer_length']) + 1):
                if len(primer_sets[prefix][k]) > 0:
                    tasks.append((primer_sets[prefix][k], k, fg_genome, data['fg_prefixes'][i]))
    immut_list = pool.map(get_positions, tasks)
    tb = pc()
    if data['verbose']:
        print("Time to get all positions:", str(round(tb - ta, 4)))

    df, primers_with_positions, chr_lens = make_df(immut_list, primer_dict, data['fg_prefixes'])

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

    return df, primers_with_positions, chr_lens

if __name__ == "__main__":
    # in_json = sys.argv[1]
    in_json = '/Users/Kaleb/Desktop/Bailey_Lab/code/newswga/params/new_params.json'
    with open(in_json, 'r') as f:
        data = json.load(f)
    df, primers_with_positions, chr_lens = main(data)