import pandas as pd
import numpy as np
import sys
import json
from numba import njit

def smith_waterman(seq1, seq2):
    match = 1
    mismatch = -1
    gap = -10
    # Generating the empty matrices for storing scores and tracing
    row = len(seq1) + 1
    col = len(seq2) + 1
    matrix = np.zeros(shape=(row, col), dtype=int)  
    
    # Initialising the variables to find the highest scoring cell
    max_score = -1
    
    # Calculating the scores for all cells in the matrix
    for i in range(1, row):
        for j in range(1, col):
            # Calculating the diagonal score (match score)
            match_value = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diagonal_score = matrix[i - 1, j - 1] + match_value
            
            # Calculating the vertical gap score
            vertical_score = matrix[i - 1, j] + gap
            
            # Calculating the horizontal gap score
            horizontal_score = matrix[i, j - 1] + gap
            
            # Taking the highest score 
            matrix[i, j] = max(0, diagonal_score, vertical_score, horizontal_score)
            
            # Tracking the cell with the maximum score
            if matrix[i, j] > max_score:
                max_score = matrix[i, j]
    return max_score

def binding_score(seq1, seq2):
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
    shorter_seq = min(len(seq1),len(seq2))
    for i in range(shorter_seq):
        score = 0
        for j in range(shorter_seq-i):
            score += comps[seq1[i+j]][seq2[j]]
        if score > max_score:
            max_score = score
    return max_score

def is_dimer(primers:list, primer_to_check:str) -> bool:
    max_alignment = 3
    reversed_primer = primer_to_check[::-1]
    if binding_score(primer_to_check, reversed_primer) > max_alignment:
        return True
    if primers == []:
        return False
    for primer in primers:
        if binding_score(primer, reversed_primer) > max_alignment:
            return True
    return False

def main(df:pd.DataFrame, primers_with_positions:dict, data):
    print("Finding sets...")
    df = df.reset_index(drop=True)
    fg_length = sum(data["fg_seq_lengths"])
    frag_length = data["fragment_length"]
    all_indices = set()
    coverage_change = 0
    primes = []
    coverage = 0.5
    index = 0
    total_fgs = 0
    total_bgs = 0
    while coverage < data["target_coverage"] and index < len(df):
        primer = df['primer'][index]
        count = df['fg_count'][index]
        if not is_dimer(primes, primer):
            old_index_count = len(all_indices)
            for pos in primers_with_positions[primer]:
                all_indices |= set(range(int(pos), min(int(pos) + frag_length, fg_length)))
            if (len(all_indices) - old_index_count) > coverage_change * count * frag_length:
                primes.append(primer)
                total_fgs += count
                total_bgs += df['bg_count'][index]
                print("Current set: " + str(primes))
                print("Current coverage: " + str((int(100000*len(all_indices))/fg_length)/100000))
        coverage = (len(all_indices))/fg_length
        index += 1
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