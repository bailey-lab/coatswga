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
                    # final check if background frequency is below threshold, add to the primer dict if yes
                    if bg_count/bg_total_length < max_freq:
                        primer_dict[kmer] = [fg_count, bg_count]
    return (prefix, primer_dict)

def make_tasks(mini, maxi, genomes, primer_dict, prefixes):
    '''
    Initalizes the tasks array so the position calculation process can be multithreaded

    Args:
        mini: Lower bound on primer length
        maxi: Upper bound on primer length
        genomes: List of foreground genomes to search
        primer_dict: Dictionary containing all of the primers that passed the filters mapping to foreground and background counts
        prefixes: List of prefixes corresponding to each foreground genome

    Returns:
        tasks: an array of immutables where each one contains a set of primers of each length in the range, a set of reverse complements of those primers, genome paths, and sequence lengths for each
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
                tasks.append((primers_per_k, reverses_per_k, k, fg_genome, prefixes[i]))
    return tasks

def get_positions(task):
    '''
    Reads through the foreground genome and records the indices at which each primer and its reverse complement occurs

    Args:
        task:
            primer_set: Set containing all primers of length k
            rev_set: Set containing all reverse complements of primer_set
            k: Length of the primers given
            fg_genome: File path of the genome to read
            prefix: Prefix of the foreground genome, used to ID the immutable 

    Returns:
        immutable:
            prefix: Prefix of the foreground genome, used to ID the immutable 
            primer_dict: Dictionary of primers of length k mapping to a list of positions where that primer is found in each chromosome for the given genome
            rev_dict: Dictionary of reverse complements of the primers mapping to a list of positions where they can be found
            lens: Dicitonary mapping chromosome numbers to their lengths
    '''
    # breaks up the task
    primer_set, rev_set, k, fg_genome, prefix = task

    # initializes the dictionaries to store the primers mapped to their positions along with the dictionary to hold the chromosome lengths
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
            # if line is chromosome header, increment chromosome number, record length of previous chromosome, intialize new dictionary to fill for that chromosome
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
            # combines lines without whitespace
            stripped = line.strip()
            combined = ''.join([prev[len(prev)-k+1:], stripped]).upper()
            # iterates through each base checking for primers and reverse complements
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
    # adds length of last chromosome 
    lens[chr] = index + k
    return (prefix, primer_dict, rev_dict, lens)

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
    revs_with_positions = {}
    chr_lens = {}
    for prefix in prefixes:
        primers_with_positions[prefix] = {}
        revs_with_positions[prefix] = {}
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
        # iterates through each chromosome in the dictionary of primers
        for chr in trip[1]:
            # adds the length of the chromosome to the above dictionary
            chr_lens[prefix][chr] = trip[3][chr]
            # if the chromosome has not been added already, initialize an empty dictionary to fill with primers
            if chr not in primers_with_positions[prefix]:
                primers_with_positions[prefix][chr] = {}
            # for each primer in the chromosome add the positions to the central dictionary
            for primer in trip[1][chr]:
                primers_with_positions[trip[0]][chr][primer] = trip[1][chr][primer]
        # iterates through each chromosome in the dictionary of reverse complements
        for chr in trip[2]:
            # if the chromosome has not been added already, initialize an empty dictionary to fill with primers
            if chr not in revs_with_positions[prefix]:
                revs_with_positions[prefix][chr] = {}
            # for each reverse complement in the chromosome add the positions to the central dictionary
            for rev in trip[2][chr]:
                revs_with_positions[trip[0]][chr][rev] = trip[2][chr][rev]
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

    tasks = make_tasks(data["min_primer_length"], data["max_primer_length"], data["fg_genomes"], primer_dict, data["fg_prefixes"])
    pool = multiprocessing.Pool(processes=int(data['cpus']))
    immut_list = pool.map(get_positions, tasks)

    df, primers_with_positions, revs_with_positions, chr_lens = make_df(immut_list, primer_dict, data['fg_prefixes'])

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
        