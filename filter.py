import json
import subprocess
import pandas as pd
import multiprocessing
import os
import melting
import numpy as np
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
    kmer_dir = data["data_dir"] + "kmer_files/"
    bg_prefix = data['bg_prefix']
    primer_set = set()
    primer_dict = {}
    subprocess.run(["kmc_tools", "-t1", "-hp", "transform", f"{kmer_dir}{prefix}_{k}mers", "reduce", f"{kmer_dir}red_{prefix}_{k}mers", f"-ci{data['min_fg_count']}"])
    subprocess.run(["kmc_tools", 
                    "-t1", 
                    "-hp", 
                    "simple", 
                    f"{kmer_dir}red_{prefix}_{k}mers", 
                    f"{kmer_dir}{bg_prefix}_{k}mers", 
                    "intersect", 
                    f"{kmer_dir}{bg_prefix}_{k}mer_counts", 
                    "-ocright"])

    # while not os.path.exists(f"{kmer_dir}{bg_prefix}_{k}mers.txt"):
    subprocess.run(["kmc_tools", "-t1", "-hp", "transform", f"{kmer_dir}{bg_prefix}_{k}mer_counts", "dump", f"{kmer_dir}{bg_prefix}_{k}mers.txt"])
    # while not os.path.exists(f"{kmer_dir}{prefix}_{k}mers.txt"):
    subprocess.run(["kmc_tools", "-t1", "-hp", "transform", f"{kmer_dir}red_{prefix}_{k}mers", "dump", f"{kmer_dir}{prefix}_{k}mers.txt"])
    subprocess.run(["rm", 
                    f"{kmer_dir}{bg_prefix}_{k}mer_counts.kmc_pre", 
                    f"{kmer_dir}{bg_prefix}_{k}mer_counts.kmc_suf", 
                    f"{kmer_dir}red_{prefix}_{k}mers.kmc_pre", 
                    f"{kmer_dir}red_{prefix}_{k}mers.kmc_suf"])
    
    with open(f"{kmer_dir}{prefix}_{k}mers.txt", 'r') as fg_file, open(f"{kmer_dir}{bg_prefix}_{k}mers.txt", 'r') as bg_file:
        kmer_to_bg = {}
        for line in bg_file:
            kmer = line.strip().split("\t")[0]
            bg_count = int(line.strip().split("\t")[1])
            kmer_to_bg[kmer] = bg_count
        for line in fg_file:
            kmer = line.strip().split("\t")[0]
            fg_count = int(line.strip().split("\t")[1])
            if kmer not in kmer_to_bg:
                kmer_to_bg[kmer] = 0
            if kmer_to_bg[kmer]/fg_count <= data['max_ratio'] and melting.temp(kmer) <= data['max_tm']:
                primer_dict[kmer] = [fg_count, kmer_to_bg[kmer]]
                primer_set.add(kmer)
                primer_set.add(rc(kmer))
    if data['verbose']:
        print(f"{len(primer_dict)} {k}-mers in {prefix}")
    if not data['write']:
        subprocess.run(["rm", f"{kmer_dir}{bg_prefix}_{k}mers.txt", f"{kmer_dir}{prefix}_{k}mers.txt"])
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
    # df['sort_val'] = df.apply(lambda x: x['fg_count']/(x['ratio'] + 0.000000001), axis=1)
    # sorted_df = df.sort_values(by='sort_val', ascending=False)

    return df, primers_with_positions, chr_lens

def bedtooler(pos_inters:dict, all_inters:dict, data_dir:str):
    '''
    Helper method to take in two lists of intervals of covered indices and use bedtools to combine them

    Args:
        pos_inters: Dictionary mapping each chromosome to a list of interval tuples of the primer to check
        all_inters: Dictionary mapping each chromosome to a list of interval tuples of the indices covered by the current primer set, or an empty 
                    dictionary so the intervals of just the current primer can be calculated
        data_dir: Path to store the intermediate pos.bed file 

    Returns:
        length: The number of indices covered by the combined ranges
        inters: List of tuples representing the combined ranges
    '''

    # initialize dictionary to hold all the intervals
    both = {}
    # iterates through each chromosome in the passed dictionary
    if len(pos_inters) >= len(all_inters):
        big = pos_inters
        small = all_inters
    else:
        big = all_inters
        small = pos_inters
    empty_counter = 0
    for chr in big:
        # if the chromosome is already in the covered indices
        if chr in small:
            both[chr] = big[chr] + small[chr]
        else:
            both[chr] = big[chr]
        # sort the list of tuples by the first value so bedtools can use it
        both[chr].sort(key=lambda x: x[0])
        if len(both[chr]) == 0:
            empty_counter += 1
    # if empty return early
    if empty_counter == len(both):
        return 0, {}
    # open the intermediate file pos[pid].bed
    pid = os.getpid()
    with open(f"{data_dir}pos{pid}.bed", 'w') as f:
        # iterate through each chromosome in the combined dictionary, writing each range to the file
        for chr in both:
            for tup in both[chr]:
                f.write(chr + "\t" + str(tup[0]) + "\t" + str(tup[1]) + "\n")

    # variables to store output
    length = 0
    inters = {}

    # outer = subprocess.run(['bedtools', 'merge', '-i', data_dir + "/pos.bed"], capture_output=True)
    # get output from bedtools, decode, strip, and split by line
    out = subprocess.check_output(['bedtools', 'merge', '-i', f"{data_dir}pos{pid}.bed"]).decode().strip().split('\n')
    # iterates through each line
    for line in out:
        # bedtools output is tab deliminated, so split on each tab
        parts = line.split('\t')
        chr = parts[0]
        # if chromosome not contained in final dictionary, intialize new list to be added to
        if chr not in inters:
            inters[chr] = []
        # start and end points of range
        start = int(parts[1])
        end = int(parts[2])
        inters[chr].append((start,end))
        # add difference to length covered
        length += end - start
    return length, inters

def cov_computer(task) -> dict:
    df, prim_pos, chr_lens, data = task
    prefixes = data['fg_prefixes']
    lens = []
    for primer in df['primer']:
        cov = 0
        for prefix in prefixes:
            prim_dict = {}
            for chr in prim_pos[prefix]:
                prim_dict[chr] = []
                if primer in prim_pos[prefix][chr]:
                    prim_dict[chr] = [(int(pos), min(chr_lens[prefix][chr], int(pos) + data['fragment_length'])) for pos in prim_pos[prefix][chr][primer]]
            length = bedtooler(prim_dict, {}, data['data_dir'])[0]
            cov += length
        lens.append(cov)
    df['cov_len'] = lens
    return df

def main(data):
    t0 = pc()

    print("Filtering primers...")
    pool = multiprocessing.Pool(processes=int(data['cpus']))
    if data['cpus'] > 1:
        tasks = []
        for prefix in data["fg_prefixes"]:
            for k in range(int(data["min_primer_length"]), int(data["max_primer_length"]) + 1):
                tasks.append((prefix, k, data))
        primer_dicts_list = pool.map(filter_primers_into_dict, tasks)
    else:
        primer_dicts_list = []
        for prefix in data["fg_prefixes"]:
            for k in range(int(data["min_primer_length"]), int(data["max_primer_length"]) + 1):
                primer_dicts_list.append(filter_primers_into_dict((prefix, k, data)))
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

    print("Calculating coverage intervals...")
    tt = pc()
    parts = np.array_split(df, data['cpus'])
    tasks = [(arr, primers_with_positions, chr_lens, data) for arr in parts]
    covers = pool.map(cov_computer, tasks)
    df_w_lens: pd.DataFrame = pd.concat([dat for dat in covers])
    df_w_lens['sort_val'] = df_w_lens.apply(lambda x: x['cov_len']/(pow(x['ratio'] + 0.000000001, 2) * x['bg_count']), axis=1)
    df_w_lens = df_w_lens.sort_values(by='sort_val', ascending=False)
    df_w_lens = df_w_lens.reset_index(drop=True)
    if data['verbose']:
        print(f"Time getting coverage lengths: {pc() - tt}")

    if data["write"]:
        df_w_lens.to_csv(data['data_dir'] + "/primers_df.csv")
        # for i, prefix in enumerate(data['fg_prefixes']):
        #     name = prefix.split('/')[-1]
        #     if len(primers_with_positions[prefix]) > 0:
        #         with open(data['data_dir'] + "/" + name + "_primers_with_positions.csv", 'w+') as f:
        #             for chr in primers_with_positions[prefix]:
        #                 f.write(">" + str(chr) + "\n")
        #                 for primer in primers_with_positions[prefix][chr]:
        #                     to_write = str(primer) + ":" + str(primers_with_positions[prefix][chr][primer]) + "\n"
        #                     f.write(to_write)

    return df_w_lens, primers_with_positions, chr_lens

if __name__ == "__main__":
    # in_json = sys.argv[1]
    in_json = '/Users/Kaleb/Desktop/Bailey_Lab/code/newswga/params/new_params.json'
    with open(in_json, 'r') as f:
        data = json.load(f)
    df, primers_with_positions, chr_lens = main(data)