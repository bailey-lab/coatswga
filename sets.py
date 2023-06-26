import pandas as pd
import numpy as np
import sys
import json
from numba import njit
from filter import rc

def binding_score(seq1, seq2):
    '''
    Helper method to calculate the "binding score" between two primers.

    Args:
        seq1: First primer to compare
        seq2: Second primer to compare

    Returns:
        max_score: The highest binding score among all possible alignments of the two primers
    '''

    # Intialized dictionary with scores for binding. +1 if the two bases bind, -1 otherwise.
    comps = {
        "A": {
            "A": -1,
            "G": -1,
            "C": -1,
            "T": 1
        },
        "G": {
            "A": -1,
            "G": -1,
            "C": 1,
            "T": -1
        },
        "C": {
            "A": -1,
            "G": 1,
            "C": -1,
            "T": -1
        },
        "T": {
            "A": 1,
            "G": -1,
            "C": -1,
            "T": -1
        }
    }
    max_score = 0

    # The length of the shorter primer is needed so there are no indexing errors
    shorter_seq = min(len(seq1),len(seq2))

    # Iterates through each possible alignment of the two primers and scores it based on the comps dictionary
    for i in range(shorter_seq):
        score = 0
        for j in range(shorter_seq-i):
            score += comps[seq1[i+j]][seq2[j]]
        if score > max_score:
            max_score = score
    return max_score

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
    max_alignment = 3 
    reversed_primer = primer_to_check[::-1]
    if primers == []:
        return False
    if binding_score(primer_to_check, reversed_primer) > max_alignment:
        return True
    for primer in primers:
        if primer_to_check in primer:
            return True
        if primer in primer_to_check:
            return True
        if binding_score(primer, reversed_primer) > max_alignment:
            return True
    return False

def main(df:pd.DataFrame, primers_with_positions:dict, data, reverse_positions:dict):
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

    print("Finding sets...")

    # indices of rows get scrambled when sorting by ratio, have to reset indices to iterate through in order
    df = df.reset_index(drop=True)

    fg_length = sum(data["fg_seq_lengths"])
    frag_length = data["fragment_length"]

    # set to hold the covered indices
    all_indices = set()

    # edits the threshold that each primer must cover a certain percent of new indices
    coverage_change = 0.9

    primes = []
    rev_index = []
    coverage = 0
    index = 0
    total_fgs = 0
    total_bgs = 0
    while coverage < data["target_coverage"] and index < len(df) and coverage_change >= 0.1:
        # Primer to check and the counts of foreground hits
        primer = df['primer'][index]

        # checks if primer forms a self dimer, is a substring of a primer already in the set, or forms a 
        # primer-primer dimer with any primers already in the set
        if primer not in primes:
            if not is_dimer(primes, primer):

                # stores the size of set before adding the next primer
                old_indices = all_indices.copy()
                new_covered = 0

                mat = [ list(filter(lambda x: x > 0 and x <= frag_length, [b - a for b in (rev_index + reverse_positions[rc(primer)])])) for a in primers_with_positions[primer]]
                primer_positions = primers_with_positions[primer]
                for i in range(len(mat)):
                    if mat[i] != []:
                        all_indices |= set(range(int(primer_positions[i]), primer_positions[i] + max(mat[i])))
                        new_covered += max(mat[i])

                # iterates through the list of positions (binding sites) of the primer, adds the covered indices to the set
                # for pos in primers_with_positions[primer]:
                #     all_indices |= set(range(int(pos), min(int(pos) + frag_length, fg_length)))

                # checks if the primer covered a high enough percent of bases that were not previously covered by the set
                if (len(all_indices) - len(old_indices)) > coverage_change * new_covered:
                    primes.append(primer)
                    rev_index = rev_index + reverse_positions[rc(primer)]
                    total_fgs += df['fg_count'][index]
                    total_bgs += df['bg_count'][index]
                    print("Current set: " + str(primes))
                    print("Current coverage: " + str((int(100000*len(all_indices))/fg_length)/100000))
                else:
                    all_indices = old_indices
        coverage = (len(all_indices))/fg_length
        index += 1
        if index == len(df):
            index = 0
            coverage_change -= 0.1
    print("\nFinal primers: " + str(primes))
    print("Expected coverage: " + str(coverage))
    print("Total foreground hits: " + str(total_fgs))
    print("Total background hits: " + str(total_bgs))
    print("Bg/fg ratio: " + str(int(10000*total_bgs/total_fgs)/10000))

if __name__ == "__main__":
    # in_json = sys.argv[1]
    in_json = '/Users/Kaleb/Desktop/Bailey_Lab/code/newswga/new_src/new_params.json'
    with open(in_json, 'r') as f:
        data = json.load(f)
    df = pd.read_csv(data['data_dir'] + "/primers_df.csv")
    prim_w_pos = {}
    with open(data['data_dir'] + "/primers_with_ginis.csv", 'r') as f:
        pair = f.readline().strip()
        while pair != "":
            primer = pair.split(':')[0]
            positions = pair.split(':')[1].strip('][').split(', ')
            prim_w_pos[primer] = positions
            pair = f.readline().strip()
    main(df, prim_w_pos, data)