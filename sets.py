import json
import pandas as pd
import multiprocessing
import os
from multiply_align.algorithms import PrimerDimerLike
from filter import rc, bedtooler
from time import perf_counter as pc

def is_dimer(primers:list, primer_to_check:str) -> bool:
    '''
    Checks if the given primer forms a self dimer, a primer-primer dimer with any other primers in the list, or is a substring
    of another primer in the list. Self dimer/primer-primer dimer score is calculated using an algorithm implemented in Multiply, developed by 
    Cesare et al. (2023), first described in Johnston et al. (2019) Sci Reports

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
        if model.score < -2.79: # threshold defined by Johnston et al. that guarantees no dimers forming
            return True
    return False 

def setter(task):
    """
    Finds a set of primers with the highest specificity that theoretically meets the target_coverage value specified in the JSON file. 
    Primers in the set should not form self dimers or primer-primer dimers within the set. Does this by treating each base in the
    foreground genome as an index, then iterates through each position of each primer in the fg genome, adding a range of indices
    (determined by the fragment_length in the JSON file) to a set tracking the covered indices. Ensures tiling by checking that the 
    percent of new indices covered surpasses a threshold determined by coverage_change.

    Args:
        task:
            primer: the primer to start with
            index: the index of the first primer in the dataframe
            df: Pandas DataFrame with all the primers that passed the filter, their foreground hits, background hits, and bg/fg ratio. The dataframe 
                is sorted first on bg/fg ratio (low to high to get best specificity) then by fg count (high to low to get highest frequency).
            primers_with_postions: A dictionary mapping chromosome number to a dictionary containing all the primers that passed the filter mapped to
                                    a list of the indices of all the hits of that primer on the fg genome
            chr_lens: Dictionary mapping chromosome number to its length
            data: A JSON object containing hyperparameters
    """

    primer, i, df, primers_with_positions, chr_lens, data = task

    pid = os.getpid()

    # gets total foreground length
    total_fg_length = 0
    for prefix in chr_lens:
        for chr in chr_lens[prefix]:
            total_fg_length += chr_lens[prefix][chr]

    # assigns variables to the fragment length and list of fg prefixes for ease of access
    frag_length = data["fragment_length"]
    prefixes = data['fg_prefixes']

    # dictionaries to hold intervals of covered indices in the forward and reverse directions
    fwd_inters = {}
    rev_inters = {}

    # dictionary to hold intervals of covered indices for each primer, speeds up program so the intervals don't have to be
    # recalculated every time the primer is checked
    prim_inters = {prefix:{} for prefix in prefixes}

    # initialzie lengths of covered indices in each direction
    fwd_len = 0
    rev_len = 0
    total_fgs = 0
    total_bgs = 0

    # edits the threshold that each primer must cover a certain percent of new indices
    cov_ch = [0.5, 0.25, 0.15, 0.1, 0.05, 0]
    coverage_change = cov_ch[0]

    primes = data["existing_primers"] + [primer]

    for prefix in prefixes:
        prim_dict = {}
        for chr in primers_with_positions[prefix]:
            if primer in primers_with_positions[prefix][chr]:
                prim_dict[chr] = [(int(pos), min(chr_lens[prefix][chr], int(pos) + frag_length)) for pos in primers_with_positions[prefix][chr][primer]]
        prim_inters[prefix][primer] = bedtooler(prim_dict, {}, data['data_dir'])
        fwd_len += prim_inters[prefix][primer][0]
        fwd_inters[prefix] = prim_inters[prefix][primer][1]

    rev = rc(primer)
    for prefix in prefixes:
        prim_dict = {}
        for chr in primers_with_positions[prefix]:
            if rev in primers_with_positions[prefix][chr]:
                prim_dict[chr] = [(max(0, int(pos) + len(rev) - frag_length), int(pos) + len(rev)) for pos in primers_with_positions[prefix][chr][rev]]
        length, inters = bedtooler(prim_dict, {}, data['data_dir'])
        rev_len += length
        rev_inters[prefix] = inters

    fwd_coverage = fwd_len / total_fg_length
    total_fgs = df['fg_count'][i]
    total_bgs = df['bg_count'][i]
    index = 0
    while fwd_coverage < data["target_coverage"] and index < len(df):
        # Primer to check and the counts of foreground hits
        primer = df['primer'][index]
        count = df['fg_count'][index]

        # checks if primer forms a self dimer, is a substring of a primer already in the set, or forms a 
        # primer-primer dimer with any primers already in the set
        if primer not in primes:
            if not is_dimer(primes, primer):
                
                # initialize variable to hold the old covered genome intervals, amount of indices covered by the primer, and amount of indices covered by the primer plus the set
                tot_len = 0
                new_inters = fwd_inters.copy()

                # for each fg genome passed
                for prefix in prefixes:
                    # check if the primer intervals have already been calculated
                    if primer not in prim_inters[prefix]:
                        # if no, run bedtooler() for the primer's intervals and an empty dictionary as the second argument, store length and intervals
                        prim_inters[prefix][primer] = {}
                        for chr in primers_with_positions[prefix]:
                            if primer in primers_with_positions[prefix][chr]:
                                prim_inters[prefix][primer][chr] = [(int(pos), min(chr_lens[prefix][chr], int(pos) + frag_length)) for pos in primers_with_positions[prefix][chr][primer]]
                    # run bedtooler for the primer and the intervals of the current set
                    length, intervals = bedtooler(prim_inters[prefix][primer], new_inters[prefix], data['data_dir'])
                    tot_len += length
                    new_inters[prefix] = intervals

                # checks if the primer covered a high enough percent of bases that were not previously covered by the set
                if (tot_len - fwd_len)/total_fg_length >= coverage_change * (1 - fwd_coverage):
                    primes.append(primer)
                    fwd_len = tot_len
                    fwd_inters = new_inters.copy()
                    coverage_change = cov_ch[0]
                    # calculate the reverse strand coverage provided by any reverse complements of the primer
                    rev_len = 0
                    rev = rc(primer)
                    for prefix in prefixes:
                        prim_dict = {}
                        for chr in primers_with_positions[prefix]:
                            if rev in primers_with_positions[prefix][chr]:
                                prim_dict[chr] = [(max(0, int(pos) + len(rev) - frag_length), int(pos) + len(rev)) for pos in primers_with_positions[prefix][chr][rev]]
                        length, inters = bedtooler(prim_dict, rev_inters[prefix], data['data_dir'])
                        rev_len += length
                        rev_inters[prefix] = inters
                    total_fgs += count
                    total_bgs += df['bg_count'][index]
                    fwd_coverage = fwd_len/total_fg_length
        index += 1
        if index == len(df):
            index = 0
            coverage_change = cov_ch[cov_ch.index(coverage_change) + 1]
            if coverage_change == 0:
                break
    
    rev_coverage = rev_len/total_fg_length
    coverage_change = cov_ch[0]
    index = 0
    prim_inters = {prefix: {} for prefix in prefixes}
    while rev_coverage < data["target_coverage"] and index < len(df) and coverage_change >= 0.01:
        # Primer to check and the counts of foreground hits
        primer = df['primer'][index]
        count = df['fg_count'][index]
        rev = rc(primer)

        # checks if primer forms a self dimer, is a substring of a primer already in the set, or forms a 
        # primer-primer dimer with any primers already in the set
        if rev not in primes:
            if not is_dimer(primes, rev):

                tot_len = 0
                new_inters = rev_inters.copy()
                for prefix in prefixes:
                    if rev not in prim_inters[prefix]:
                        prim_inters[prefix][rev] = {}
                        for chr in primers_with_positions[prefix]:
                            if primer in primers_with_positions[prefix][chr]:
                                prim_inters[prefix][rev][chr] = [(max(0, int(pos) + len(primer) - frag_length), int(pos) + len(primer)) for pos in primers_with_positions[prefix][chr][primer]]
                    length, intervals = bedtooler(prim_inters[prefix][rev], new_inters[prefix], data['data_dir'])
                    tot_len += length
                    new_inters[prefix] = intervals 

                # checks if the primer covered a high enough percent of bases that were not previously covered by the set
                if (tot_len - rev_len)/total_fg_length >= coverage_change * (1 - rev_coverage):
                    primes.append(rev)
                    rev_len = tot_len
                    rev_inters = new_inters.copy()
                    coverage_change = cov_ch[0]
                    if fwd_coverage < data["target_coverage"]:
                        fwd_len = 0
                        for prefix in prefixes:
                            prim_dict = {}
                            for chr in primers_with_positions[prefix]:
                                if primer in primers_with_positions[prefix][chr]:
                                    prim_dict[chr] = [(int(pos), min(chr_lens[prefix][chr], int(pos) + frag_length)) for pos in primers_with_positions[prefix][chr][primer]]
                            length, intervals = bedtooler(prim_dict, new_inters[prefix], data['data_dir'])
                            fwd_len += length
                            fwd_inters[prefix] = intervals
                    fwd_coverage = fwd_len/total_fg_length
                    total_fgs += count
                    total_bgs += df['bg_count'][index]
                    rev_coverage = rev_len/total_fg_length
        index += 1
        if index == len(df):
            index = 0
            coverage_change = cov_ch[cov_ch.index(coverage_change) + 1]
            if coverage_change == 0:
                break
    
    index = 0
    while len(primes) < data['min_set_size'] and index < len(df):
        primer = df['primer'][index]
        if not is_dimer(primes, primer):
            primes.append(primer)
            total_fgs += df['fg_count'][index]
            total_bgs += df['bg_count'][index]
        index += 1
    
    if data['verbose']:
        print(f"\nFinal set for {pid}: ")
        print(str(primes))
        print("Expected forward coverage: " + str(round(fwd_coverage, 3)))
        print("Expected reverse coverage: " + str(round(rev_coverage, 3)))
        print("Total foreground hits: " + str(total_fgs))
        print("Total background hits: " + str(total_bgs))
        print("Bg/fg ratio: " + str(round(total_bgs/total_fgs, 3)) + "\n")

    return (primes, fwd_coverage, rev_coverage, total_fgs, total_bgs, total_bgs/total_fgs)

def main(df:list, primers_with_positions:dict, chr_lens:dict, data):
    print("Finding sets...")
    t0 = pc()
    pool = multiprocessing.Pool(processes=data['cpus'])

    if data['force_coverage_threshold']:
        thresh = data['target_coverage']
    else:
        thresh = 0.01
    max_cov = 0
    index = 0
    all_out = []
    while max_cov < thresh and (index + data['cpus']) < len(df):
        tasks = []
        for i in range(index, index + data['cpus']):
            if not is_dimer(data["existing_primers"], df['primer'][i]):
                tasks.append((df['primer'][i], i, df, primers_with_positions, chr_lens, data))
        out = pool.map(setter, tasks)
        all_out.extend(out)
        index += data['cpus']
        for outer in out:
            av = (outer[1] + outer[2])/2
            if av > max_cov:
                max_cov = av

    best_coverage = 0
    coverage = (all_out[0][1] + all_out[0][2])/2
    lowest_ratio = 0
    ratio = all_out[0][5]
    for i in range(1, len(all_out)):
        if (all_out[i][1] + all_out[i][2])/2 > coverage:
            best_coverage = i
            coverage = (all_out[i][1] + all_out[i][2])/2
        elif (all_out[i][1] + all_out[i][2])/2 == coverage and all_out[i][5] < all_out[best_coverage][5]:
            best_coverage = i
        
        if all_out[i][5] < ratio:
            lowest_ratio = i
            ratio = all_out[i][5]
        elif all_out[i][5] == ratio and (all_out[i][1] + all_out[i][2])/2 > (all_out[lowest_ratio][1] + all_out[lowest_ratio][2])/2:
            lowest_ratio = i
    
    if best_coverage == lowest_ratio:
        print("\nSet with highest coverage, fewest primers, and lowest ratio: ")
        print(str(all_out[best_coverage][0]))
        print("Expected forward coverage: " + str(round(all_out[best_coverage][1], 3)))
        print("Expected reverse coverage: " + str(round(all_out[best_coverage][2], 3)))
        print("Total foreground hits: " + str(all_out[best_coverage][3]))
        print("Total background hits: " + str(all_out[best_coverage][4]))
        print("Bg/fg ratio: " + str(round(all_out[best_coverage][5], 3)))
    else:
        print("\nSet with lowest ratio: ")
        print(str(all_out[lowest_ratio][0]))
        print("Expected forward coverage: " + str(round(all_out[lowest_ratio][1], 3)))
        print("Expected reverse coverage: " + str(round(all_out[lowest_ratio][2], 3)))
        print("Total foreground hits: " + str(all_out[lowest_ratio][3]))
        print("Total background hits: " + str(all_out[lowest_ratio][4]))
        print("Bg/fg ratio: " + str(round(all_out[lowest_ratio][5], 3)))

        print("\nSet with highest coverage: ")
        print(str(all_out[best_coverage][0]))
        print("Expected forward coverage: " + str(round(all_out[best_coverage][1], 3)))
        print("Expected reverse coverage: " + str(round(all_out[best_coverage][2], 3)))
        print("Total foreground hits: " + str(all_out[best_coverage][3]))
        print("Total background hits: " + str(all_out[best_coverage][4]))
        print("Bg/fg ratio: " + str(round(all_out[best_coverage][5], 3)) + "\n")
    os.system(f"rm {data['data_dir']}pos*.bed")
    print("Time finding sets:", pc() - t0)

if __name__ == "__main__":
    # in_json = sys.argv[1]
    in_json = '/Users/Kaleb/Desktop/Bailey_Lab/code/newswga/params/soaper_params.json'
    with open(in_json, 'r') as f:
        data = json.load(f)
    # df = pd.read_csv(data['data_dir'] + "/primers_df.csv")
    df = pd.read_csv("/Users/kaleb/Desktop/soap_counts.csv")
    df = df.sort_values(by=["ratio", "fg_count"], ascending=[True, False])
    pref = '/Users/Kaleb/Desktop/Bailey_Lab/code/newswga/kmers/malaria'
    prim_w_pos = {pref: {}}
    # with open(data['data_dir'] + "/primers_with_positions.csv", 'r') as f:
    #     pair = f.readline().strip()
    #     while pair != "":
    #         primer = pair.split(':')[0]
    #         positions = pair.split(':')[1].strip('][').split(', ')
    #         if positions == ['']:
    #             positions = []
    #         prim_w_pos[primer] = positions
    #         pair = f.readline().strip()
    with open("/Users/kaleb/Desktop/soap_primers_with_positions.csv", 'r') as f:
        chr = ""
        for line in f:
            line = line.strip()
            if ">" in line:
                chr = line[1:]
                prim_w_pos[pref][chr] = {}
            else:
                primer = line.split(':')[0]
                positions = line.split(':')[1].strip('][').split(', ')
                if positions == ['']:
                    positions = []
                prim_w_pos[pref][chr][primer] = positions
    lens = {pref:{}}
    with open("/Users/kaleb/Desktop/soap_chr_lens.csv", "r") as f:
        for line in f:
            line = line.strip().split("~")
            lens[pref][line[0]] = int(line[1])
    main(df, prim_w_pos, lens, data)