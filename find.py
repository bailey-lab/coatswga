import sys
import os
import json

def run_jellyfish(data, genome_fname=None, output_prefix=None, min=6, max=12):
    """
    Runs jellyfish program using the output_prefix and transfroms the kmer count information txt files. Count k-mers from 6 to 12.

    Args:
        genome_fname: The fasta file used to count kmers.
        output_prefix: The output path prefix for the output files. Resulting output files will be suffixed by _kmer_all.txt for k from 6 to 12 inclusive.
    """
    for k in range(min, max+1, 1):
        if not os.path.exists(output_prefix+'_'+str(k)+'mer_all.txt'):
            os.system("jellyfish count -m "+str(k) + " -s 1000000 -t " + str(data['cpus']) + " " + genome_fname + " -o " + output_prefix+'_'+str(k)+'mer_all.jf')
            os.system("jellyfish dump -c " + output_prefix+'_'+str(k)+'mer_all.jf' + " > " + output_prefix+'_'+str(k)+'mer_all.txt')
        if os.path.exists(output_prefix+'_'+str(k)+'mer_all.jf'):
            os.system("rm " + output_prefix+'_'+str(k)+'mer_all.jf')

def step1(data, fg_prefixes, fg_genomes, bg_prefixes, bg_genomes, min, max):
    """
    Creates files of all k-mers of specified lengths which is located in the path specificed by --kmer_fore
    """
    for prefix in fg_prefixes:
        if not os.path.exists(os.path.dirname(prefix)):
            os.makedirs(os.path.dirname(prefix))

    print("Running jellyfish for foreground...")
    for i, fg_prefix in enumerate(fg_prefixes):
        run_jellyfish(data, fg_genomes[i], fg_prefix, min, max)

    print("Running jellyfish for background...")
    for i, bg_prefix in enumerate(bg_prefixes):
        for k in range(min, max+1, 1):
            if not os.path.exists(bg_prefix+'_'+str(k)+'mer_all.txt'):
                os.system("jellyfish count -m "+str(k) + " -s 1000000 -t " + str(data['cpus']) + " " + bg_genomes[i] + " -o " + bg_prefix+'_'+str(k)+'mer_all.jf')

    print("Done running jellyfish")

def main(data):
    step1(data, data["fg_prefixes"], data["fg_genomes"], data['bg_prefixes'], data['bg_genomes'], int(data["min_primer_length"]), int(data["max_primer_length"]))

if __name__ == "__main__":
    in_json = sys.argv[1]
    with open(in_json, 'r') as f:
        global data
        data = json.load(f)
    main(data)