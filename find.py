import json
import subprocess
import os

def step1(data_dir, fg_prefixes, fg_genomes, bg_prefix, bg_genomes, min, max, cpus):
    """
    Creates files of all k-mers of specified lengths which is located in the path specificed in the JSON file
    """
    if not os.path.exists(os.path.dirname(data_dir)):
        os.makedirs(os.path.dirname(data_dir))

    kmer_dir = data_dir + "kmer_files/"
    if not os.path.exists(os.path.dirname(kmer_dir)):
        os.makedirs(os.path.dirname(kmer_dir))


    print("Running kmc for foreground...")
    for i, fg_prefix in enumerate(fg_prefixes):
        for k in range(min, max+1, 1):
            if not os.path.exists(kmer_dir + fg_prefix+'_'+str(k)+'mers.kmc_pre') or not os.path.exists(kmer_dir + fg_prefix+'_'+str(k)+'mers.kmc_suf'):
                subprocess.run(["kmc", f"-k{k}", "-hp", f"-t{cpus}", "-fm", "-cs1000000000", "-b", f"{fg_genomes[i]}", f"{kmer_dir}{fg_prefix}_{k}mers", f"{kmer_dir}"], stdout=subprocess.DEVNULL)

    print("Running kmc for background...")
    with open(kmer_dir + "/files", 'w') as f:
        for genome in bg_genomes:
            f.write(genome + "\n")
    for k in range(min, max+1, 1):
        if not os.path.exists(f'{kmer_dir}{bg_prefix}_{k}mers.kmc_pre') or not os.path.exists(f'{kmer_dir}{bg_prefix}_{k}mers.kmc_suf'):
            subprocess.run(["kmc", f"-k{k}", "-hp", f"-t{cpus}", "-fm", "-cs1000000000", "-b", f"@{kmer_dir}/files", f"{kmer_dir}{bg_prefix}_{k}mers", f"{kmer_dir}"], stdout=subprocess.DEVNULL)
    os.system(f"rm {kmer_dir}/files")

    print("Done running kmc")

def main(data):
    step1(data['data_dir'], data["fg_prefixes"], data["fg_genomes"], data['bg_prefix'], data['bg_genomes'], int(data["min_primer_length"]), int(data["max_primer_length"]), data['cpus'])

if __name__ == "__main__":
    in_json = '/Users/kaleb/Desktop/Bailey_Lab/code/newswga/params/new_params.json'
    with open(in_json, 'r') as f:
        global data
        data = json.load(f)
    main(data)