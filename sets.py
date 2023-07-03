import pandas as pd
import numpy as np
import multiprocessing as mp
import json
from numba import njit
from multiply_align.algorithms import PrimerDimerLike
from filter import rc

def get_coverages(task):
    df, prim_w_pos, frag_length, fg_length = task
    fwd_cov = {}
    rev_cov = {}
    for primer in df:
        fwd = set()
        rev = set()
        for pos in prim_w_pos[primer]:
            fwd |= set(range(int(pos), min(int(pos) + frag_length, fg_length)))
            rev |= set(range(max(int(pos) - frag_length, 0), int(pos)))
        fwd_cov[primer] = fwd
        rev_cov[primer] = rev
    return (fwd_cov, rev_cov)

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

def main(df:pd.DataFrame, primers_with_positions:dict, rev_positions, data):
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

    fwd_sets = {}
    rev_sets = {}
    print("Calculating coverage sets...")
    if data['cpus'] > 1:
        tasks = []
        for sub in np.array_split(df['primer'], data['cpus']):
            tasks.append((sub, primers_with_positions, frag_length, fg_length))
        pool = mp.Pool(processes=len(tasks))
        coverage_tups = pool.map(get_coverages, tasks)
        for tup in coverage_tups:
            fwd_sets.update(tup[0])
            rev_sets.update(tup[1])
    else:
        fwd_sets, rev_sets = get_coverages(df['primer'], primers_with_positions, data)
    


    # set to hold the covered indices
    fwd_indices = set()
    rev_indices = set()

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

                # stores the size of set before adding the next primer
                old_fwd = fwd_indices.copy()

                # iterates through the list of positions (binding sites) of the primer, adds the covered indices to the set
                fwd_indices |= fwd_sets[primer]

                # checks if the primer covered a high enough percent of bases that were not previously covered by the set
                if (len(fwd_indices) - len(old_fwd)) > coverage_change * len(fwd_sets[primer]):
                    primes.append(primer)
                    for pos in rev_positions[rc(primer)]:
                        rev_indices |= set(range(max(int(pos) - frag_length, 0), int(pos)))
                    total_fgs += count
                    total_bgs += df['bg_count'][index]
                    fwd_coverage = len(fwd_indices)/fg_length
                    print(str(primer) + " added to set")
                    print("Current forward coverage: " + str((int(1000*len(fwd_indices))/fg_length)/1000))
                else:
                    fwd_indices = old_fwd
        index += 1
        if index == len(df):
            index = 0
            coverage_change -= 0.1
    
    rev_coverage = len(rev_indices)/fg_length
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

                # stores the size of set before adding the next primer
                old_rev = rev_indices.copy()

                # iterates through the list of positions (binding sites) of the primer, adds the covered indices to the set
                rev_indices |= rev_sets[primer]

                # checks if the primer covered a high enough percent of bases that were not previously covered by the set
                if (len(rev_indices) - len(old_rev)) > coverage_change * len(rev_sets[primer]):
                    primes.append(rc(primer))
                    total_fgs += count
                    total_bgs += df['bg_count'][index]
                    rev_coverage = len(rev_indices)/fg_length
                    # print("Current set: " + str(primes))
                    # print("Current reverse coverage: " + str((int(1000*len(rev_indices))/fg_length)/1000))
                else:
                    rev_indices = old_rev
        index += 1
        if index == len(df):
            index = 0
            coverage_change -= 0.1

    print("\nFinal primers: " + str(primes))
    print("Expected forward coverage: " + str(fwd_coverage))
    print("Expected reverse coverage: " + str(rev_coverage))
    print("Total foreground hits: " + str(total_fgs))
    print("Total background hits: " + str(total_bgs))
    print("Bg/fg ratio: " + str(int(1000*total_bgs/total_fgs)/1000))

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