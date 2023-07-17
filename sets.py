import pandas as pd
import json
from multiply_align.algorithms import PrimerDimerLike
from filter import rc
import subprocess
import multiprocessing
import os
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
    for chr in pos_inters:
        # if the chromosome is already in the covered indices
        if chr in all_inters:
            both[chr] = pos_inters[chr] + all_inters[chr]
        else:
            both[chr] = pos_inters[chr]
        # sort the list of tuples by the first value so bedtools can use it
        both[chr].sort(key=lambda x: x[0])
    # if empty return early
    if both == {}:
        return 0, {}
    # open the intermediate file pos[pid].bed
    pid = os.getpid()
    with open(f"{data_dir}/pos{pid}.bed", 'w') as f:
        # iterate through each chromosome in the combined dictionary, writing each range to the file
        for chr in both:
            for tup in both[chr]:
                f.write(chr + "\t" + str(tup[0]) + "\t" + str(tup[1]) + "\n")

    # variables to store output
    length = 0
    inters = {}

    # outer = subprocess.run(['bedtools', 'merge', '-i', data_dir + "/pos.bed"], capture_output=True)
    # get output from bedtools, decode, strip, and split by line
    out = subprocess.check_output(['bedtools', 'merge', '-i', f"{data_dir}/pos{pid}.bed"]).decode().strip().split('\n')
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

    primer, index, df, primers_with_positions, chr_lens, data = task

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

    # edits the threshold that each primer must cover a certain percent of new indices
    coverage_change = 0.9

    primes = [primer]

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
    total_fgs = df['fg_count'][index]
    total_bgs = df['bg_count'][index]
    while fwd_coverage < data["target_coverage"] and index < len(df) and coverage_change >= 0.1:
        # Primer to check and the counts of foreground hits
        primer = df['primer'][index]
        count = df['fg_count'][index]

        # checks if primer forms a self dimer, is a substring of a primer already in the set, or forms a 
        # primer-primer dimer with any primers already in the set
        if primer not in primes:
            if not is_dimer(primes, primer):
                
                # initialize variable to hold the old covered genome intervals, amount of indices covered by the primer, and amount of indices covered by the primer plus the set
                tot_len = 0
                prim_coverage_len = 0
                new_inters = fwd_inters.copy()

                # for each fg genome passed
                for prefix in prefixes:
                    # check if the primer intervals have already been calculated
                    if primer not in prim_inters[prefix]:
                        # if no, run bedtooler() for the primer's intervals and an empty dictionary as the second argument, store length and intervals
                        prim_dict = {}
                        for chr in primers_with_positions[prefix]:
                            if primer in primers_with_positions[prefix][chr]:
                                prim_dict[chr] = [(int(pos), min(chr_lens[prefix][chr], int(pos) + frag_length)) for pos in primers_with_positions[prefix][chr][primer]]
                        prim_inters[prefix][primer] = bedtooler(prim_dict, {}, data['data_dir'])
                    # run bedtooler for the primer and the intervals of the current set
                    length, intervals = bedtooler(prim_inters[prefix][primer][1], new_inters[prefix], data['data_dir'])
                    prim_coverage_len += prim_inters[prefix][primer][0]
                    tot_len += length
                    new_inters[prefix] = intervals

                # checks if the primer covered a high enough percent of bases that were not previously covered by the set
                if (tot_len - fwd_len) >= coverage_change * prim_coverage_len:
                    primes.append(primer)
                    fwd_len = tot_len
                    fwd_inters = new_inters.copy()

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
                    if data['verbose']:
                        print(f"{pid} added  {primer} to set")
                        print(f"{pid} forward coverage: " + str(round(fwd_coverage, 3)))
        index += 1
        if index == len(df):
            index = 0
            coverage_change = round(coverage_change - 0.1, 2)
            if data['verbose']:
                print(f"{pid} coverage factor to {coverage_change}")
    
    rev_coverage = rev_len/total_fg_length
    coverage_change = 0.9
    index = 0
    prim_inters = {prefix: {} for prefix in prefixes}
    while rev_coverage < data["target_coverage"] and index < len(df) and coverage_change >= 0.1:
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
                prim_coverage_len = 0
                for prefix in prefixes:
                    if rev not in prim_inters[prefix]:
                        prim_dict = {}
                        for chr in primers_with_positions[prefix]:
                            if primer in primers_with_positions[prefix][chr]:
                                prim_dict[chr] = [(max(0, int(pos) + len(primer) - frag_length), int(pos) + len(primer)) for pos in primers_with_positions[prefix][chr][primer]]
                        prim_inters[prefix][rev] = bedtooler(prim_dict, {}, data['data_dir'])
                    prim_coverage_len += prim_inters[prefix][rev][0]

                    length, intervals = bedtooler(prim_inters[prefix][rev][1], new_inters[prefix], data['data_dir'])
                    tot_len += length
                    new_inters[prefix] = intervals 

                # checks if the primer covered a high enough percent of bases that were not previously covered by the set
                if (tot_len - rev_len) >= coverage_change * prim_coverage_len:
                    primes.append(rev)
                    rev_len = tot_len
                    rev_inters = new_inters.copy()
                    total_fgs += count
                    total_bgs += df['bg_count'][index]
                    rev_coverage = rev_len/total_fg_length
                    if data['verbose']:
                        print(f"{pid} added  {primer} to set")
                        print(f"{pid} reverse coverage: " + str(round(rev_coverage, 3)))
        index += 1
        if index == len(df):
            index = 0
            coverage_change = round(coverage_change - 0.1, 2)
            if data['verbose']:
                print(f"{pid} coverage factor to {coverage_change}")
    if os.path.exists(f"{data['data_dir']}/pos{pid}.bed"):
        os.remove(f"{data['data_dir']}/pos{pid}.bed")
    
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
    # indices of rows get scrambled when sorting by ratio, have to reset indices to iterate through in order
    df = df.reset_index(drop=True)

    

    # if data['cpus'] == 1:
    #     out = setter((df['primer'][0], 0, df, primers_with_positions, chr_lens, data))
    #     print("\nFinal set: ")
    #     print(str(out[0]))
    #     print("Expected forward coverage: " + str(round(out[1], 3)))
    #     print("Expected reverse coverage: " + str(round(out[2], 3)))
    #     print("Total foreground hits: " + str(out[3]))
    #     print("Total background hits: " + str(out[4]))
    #     print("Bg/fg ratio: " + str(round(out[5], 3)) + "\n")
    # else: 
    #     tasks = []
    #     for i in range(data['cpus']):
    #         tasks.append((df['primer'][i], i, df, primers_with_positions, chr_lens, data))
    #     pool = multiprocessing.Pool(processes=data['cpus'])
    #     out = pool.map(setter, tasks)
    #     fewest = 0
    #     count = len(out[0][0])
    #     best_coverage = 0
    #     coverage = (out[0][1] + out[0][2])/2
    #     lowest_ratio = 0
    #     ratio = out[0][5]
    #     for i in range(1, len(out)):
    #         if len(out[i][0]) < count:
    #             fewest = i
    #             count = len(out[i][0])
    #         elif len(out[i][0]) == count and out[i][5] < out[fewest][5]:
    #             fewest = i
            
    #         if (out[i][1] + out[i][2])/2 > coverage:
    #             best_coverage = i
    #             coverage = (out[i][1] + out[i][2])/2
    #         elif (out[i][1] + out[i][2])/2 == coverage and out[i][5] < out[best_coverage][5]:
    #             best_coverage = i
            
    #         if out[i][5] < ratio:
    #             lowest_ratio = i
    #             ratio = out[i][5]
    #         elif out[i][5] == ratio and (out[i][1] + out[i][2])/2 > (out[lowest_ratio][1] + out[lowest_ratio][2])/2:
    #             lowest_ratio = i
    
    #     if fewest == best_coverage:
    #         print("\nSet with highest coverage and fewest primers: ")
    #         print(str(out[fewest][0]))
    #         print("Expected forward coverage: " + str(round(out[fewest][1], 3)))
    #         print("Expected reverse coverage: " + str(round(out[fewest][2], 3)))
    #         print("Total foreground hits: " + str(out[fewest][3]))
    #         print("Total background hits: " + str(out[fewest][4]))
    #         print("Bg/fg ratio: " + str(round(out[fewest][5], 3)))

    #         print("\nSet with lowest ratio: ")
    #         print(str(out[lowest_ratio][0]))
    #         print("Expected forward coverage: " + str(round(out[lowest_ratio][1], 3)))
    #         print("Expected reverse coverage: " + str(round(out[lowest_ratio][2], 3)))
    #         print("Total foreground hits: " + str(out[lowest_ratio][3]))
    #         print("Total background hits: " + str(out[lowest_ratio][4]))
    #         print("Bg/fg ratio: " + str(round(out[lowest_ratio][5], 3)) + "\n")
    #     elif fewest == lowest_ratio:
    #         print("\nSet with lowest ratio and fewest primers: ")
    #         print(str(out[fewest][0]))
    #         print("Expected forward coverage: " + str(round(out[fewest][1], 3)))
    #         print("Expected reverse coverage: " + str(round(out[fewest][2], 3)))
    #         print("Total foreground hits: " + str(out[fewest][3]))
    #         print("Total background hits: " + str(out[fewest][4]))
    #         print("Bg/fg ratio: " + str(round(out[fewest][5], 3)))

    #         print("\nSet with highest coverage: ")
    #         print(str(out[best_coverage][0]))
    #         print("Expected forward coverage: " + str(round(out[best_coverage][1], 3)))
    #         print("Expected reverse coverage: " + str(round(out[best_coverage][2], 3)))
    #         print("Total foreground hits: " + str(out[best_coverage][3]))
    #         print("Total background hits: " + str(out[best_coverage][4]))
    #         print("Bg/fg ratio: " + str(round(out[best_coverage][5], 3)) + "\n")
    #     elif best_coverage == lowest_ratio:
    #         print("\nSet with lowest ratio and highest coverage: ")
    #         print(str(out[lowest_ratio][0]))
    #         print("Expected forward coverage: " + str(round(out[lowest_ratio][1], 3)))
    #         print("Expected reverse coverage: " + str(round(out[lowest_ratio][2], 3)))
    #         print("Total foreground hits: " + str(out[lowest_ratio][3]))
    #         print("Total background hits: " + str(out[lowest_ratio][4]))
    #         print("Bg/fg ratio: " + str(round(out[lowest_ratio][5], 3)))

    #         print("\nSet with fewest primers: ")
    #         print(str(out[fewest][0]))
    #         print("Expected forward coverage: " + str(round(out[fewest][1], 3)))
    #         print("Expected reverse coverage: " + str(round(out[fewest][2], 3)))
    #         print("Total foreground hits: " + str(out[fewest][3]))
    #         print("Total background hits: " + str(out[fewest][4]))
    #         print("Bg/fg ratio: " + str(round(out[fewest][5], 3)) + "\n")
    #     else:
    #         print("\nSet with fewest primers: ")
    #         print(str(out[fewest][0]))
    #         print("Expected forward coverage: " + str(round(out[fewest][1], 3)))
    #         print("Expected reverse coverage: " + str(round(out[fewest][2], 3)))
    #         print("Total foreground hits: " + str(out[fewest][3]))
    #         print("Total background hits: " + str(out[fewest][4]))
    #         print("Bg/fg ratio: " + str(round(out[fewest][5], 3)))

    #         print("\nSet with lowest ratio: ")
    #         print(str(out[lowest_ratio][0]))
    #         print("Expected forward coverage: " + str(round(out[lowest_ratio][1], 3)))
    #         print("Expected reverse coverage: " + str(round(out[lowest_ratio][2], 3)))
    #         print("Total foreground hits: " + str(out[lowest_ratio][3]))
    #         print("Total background hits: " + str(out[lowest_ratio][4]))
    #         print("Bg/fg ratio: " + str(round(out[lowest_ratio][5], 3)))

    #         print("\nSet with highest coverage: ")
    #         print(str(out[best_coverage][0]))
    #         print("Expected forward coverage: " + str(round(out[best_coverage][1], 3)))
    #         print("Expected reverse coverage: " + str(round(out[best_coverage][2], 3)))
    #         print("Total foreground hits: " + str(out[best_coverage][3]))
    #         print("Total background hits: " + str(out[best_coverage][4]))
    #         print("Bg/fg ratio: " + str(round(out[best_coverage][5], 3)) + "\n")
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