import pandas as pd
import numpy as np
import json
from numba import njit
from multiply_align.algorithms import PrimerDimerLike
from filter import rc
import subprocess

def is_dimer(primers:list, primer_to_check:str) -> bool:
    '''
    Checks if the given primer forms a self dimer, a primer-primer dimer with any other primers in the list, or is a substring
    of another primer in the list. Self dimer/primer-primer dimer score is determined by the amount of binding bases between the 
    two.

    Args:
        primers: list of the current set of primers
        primer_to_check: primer to check for dimerization or repitition

    Returns:
        bool: True if the primer forms a self dimer, primer-primer dimer, or is already contained. False otherwise
    '''
    model = PrimerDimerLike()
    model.load_parameters()
    if primers == []:
        return False
    for primer in primers:
        if primer_to_check in primer or primer in primer_to_check:
            return True
        model.set_primers(primer, primer_to_check)
        model.align()
        if model.score < -2.79:
            return True
    return False

def bedtooler(pos_tups:list, all_inters:list, data_dir:str):
    length = 0
    inters = []

    both = pos_tups + all_inters
    both.sort(key=lambda x: x[0])
    with open(data_dir + "/pos.bed", 'w') as f:
        for tup in both:
            f.write("chr1\t" + str(tup[0]) + "\t" + str(tup[1]) + "\n")
    # outer = subprocess.run(['bedtools', 'merge', '-i', data_dir + "/pos.bed"], capture_output=True)
    out = subprocess.check_output(['bedtools', 'merge', '-i', data_dir + "/pos.bed"]).decode().strip().split('\n')
    for line in out:
        parts = line.split('\t')
        chr = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        inters.append((start,end))
        length += end - start
    return length, inters

def main(df:pd.DataFrame, primers_with_positions:dict, rev_positions:dict, data):
    """
    Finds a set of primers with the highest specificity that theoretically meets the target_coverage value specified in the JSON file. 
    Primers in the set should not form self dimers or primer-primer dimers within the set. Does this by treating each base in the
    foreground genome as an index, then iterates through each position of each primer in the fg genome, adding a range of indices
    (determined by the fragment_length in the JSON file) to a set tracking the covered indices. Ensures tiling by checking that the 
    percent of new indices covered surpasses a threshold determined by coverage_change.

    Args:
        df: Pandas DataFrame with all the primers that passed the filter, their foreground hits, background hits, gini index, and 
            bg/fg ratio. The dataframe is sorted based on bg/fg ratio, so the primers at the top have more hits on the foreground vs
            background.
        primers_with_postions: A dictionary containing all the primers that passed the filter mapped to a list of the indices of all the 
            hits of that primer on the fg genome
        data: A JSON object containing hyperparameters
    """

    # indices of rows get scrambled when sorting by ratio, have to reset indices to iterate through in order
    df = df.reset_index(drop=True)

    fg_length = sum(data["fg_seq_lengths"])
    frag_length = data["fragment_length"]

    # set to hold the covered indices
    fwd_inters = []
    fwd_len = 0
    rev_inters = []
    rev_len = 0

    # edits the threshold that each primer must cover a certain percent of new indices
    coverage_change = 0.9

    primes = []
    fwd_coverage = 0
    index = 0
    total_fgs = 0
    total_bgs = 0
    print("Finding forward set...")
    while fwd_coverage < data["target_coverage"] and index < len(df) and coverage_change >= 0.1:
        # Primer to check and the counts of foreground hits
        primer = df['primer'][index]
        count = df['fg_count'][index]

        # checks if primer forms a self dimer, is a substring of a primer already in the set, or forms a 
        # primer-primer dimer with any primers already in the set
        if primer not in primes:
            if not is_dimer(primes, primer):

                length, intervals = bedtooler([(int(pos), min(fg_length, int(pos) + frag_length)) for pos in primers_with_positions[primer]], 
                                              fwd_inters, data['data_dir'])

                # checks if the primer covered a high enough percent of bases that were not previously covered by the set
                if (length - fwd_len) > coverage_change * count * frag_length:
                    primes.append(primer)
                    fwd_len = length
                    fwd_inters = intervals
                    if rev_positions[rc(primer)] != ['']:
                        rev_len, rev_inters = bedtooler([(max(0, int(pos) - frag_length), int(pos)) for pos in rev_positions[rc(primer)]], rev_inters, data['data_dir'])
                    total_fgs += count
                    total_bgs += df['bg_count'][index]
                    fwd_coverage = fwd_len/fg_length
                    print(str(primer) + " added to set")
                    print("Current forward coverage: " + str(round(fwd_coverage, 3)))
        index += 1
        if index == len(df):
            index = 0
            coverage_change = round(coverage_change - 0.1, 2)
            print("Coverage change factor reduced to " + str(coverage_change))
    
    rev_coverage = rev_len/fg_length
    coverage_change = 0.9
    index = 0
    print("Finding reverse set...")
    while rev_coverage < data["target_coverage"] and index < len(df) and coverage_change >= 0.1:
        # Primer to check and the counts of foreground hits
        primer = df['primer'][index]
        count = df['fg_count'][index]

        # checks if primer forms a self dimer, is a substring of a primer already in the set, or forms a 
        # primer-primer dimer with any primers already in the set
        if primer not in primes:
            if not is_dimer(primes, rc(primer)):

                length, intervals = bedtooler([(max(0, int(pos) - frag_length), int(pos)) for pos in primers_with_positions[primer]], 
                                              rev_inters, data['data_dir'])
                
                # checks if the primer covered a high enough percent of bases that were not previously covered by the set
                if (length - rev_len) > coverage_change * count * frag_length:
                    primes.append(rc(primer))
                    total_fgs += count
                    total_bgs += df['bg_count'][index]
                    rev_coverage = rev_len/fg_length
                    print(str(primer) + " added to set")
                    print("Current reverse coverage: " + str(round(rev_coverage, 3)))
        index += 1
        if index == len(df):
            index = 0
            coverage_change = round(coverage_change - 0.1, 2)

    print("\nFinal primers: " + str(primes))
    print("Expected forward coverage: " + str(fwd_coverage))
    print("Expected reverse coverage: " + str(rev_coverage))
    print("Total foreground hits: " + str(total_fgs))
    print("Total background hits: " + str(total_bgs))
    print("Bg/fg ratio: " + str(round(total_bgs/total_fgs, 3)))

if __name__ == "__main__":
    # in_json = sys.argv[1]
    in_json = '/Users/Kaleb/Desktop/Bailey_Lab/code/newswga/params/new_params.json'
    with open(in_json, 'r') as f:
        data = json.load(f)
    df = pd.read_csv(data['data_dir'] + "/primers_df.csv")
    prim_w_pos = {}
    with open(data['data_dir'] + "/primers_with_positions.csv", 'r') as f:
        pair = f.readline().strip()
        while pair != "":
            primer = pair.split(':')[0]
            positions = pair.split(':')[1].strip('][').split(', ')
            prim_w_pos[primer] = positions
            pair = f.readline().strip()
    rev_w_pos = {}
    with open(data['data_dir'] + "/reverses_with_positions.csv", 'r') as f:
        pair = f.readline().strip()
        while pair != "":
            primer = pair.split(':')[0]
            positions = pair.split(':')[1].strip('][').split(', ')
            rev_w_pos[primer] = positions
            pair = f.readline().strip()
    main(df, prim_w_pos, rev_w_pos, data)