import os
import sys
import re
import subprocess
import shlex
import pandas as pd
from numpy import log10
import gzip
from random import shuffle

# Needs module load minimap2/2.24 samtools/1.9 qualimap/2.2.1 java/jdk-11 R/3.5.2
# also needs conda environment with bamcov (deeptools) and clair3

def counter(path_to_files: str, library_name: str):
    """Counts bases in a single merge_{name}.fastq.gz file."""
    count = 0
    file_path = os.path.join(path_to_files, f"{library_name}.fastq.gz")
    
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found.")
        return 0

    # Counting bases line by line to save memory
    with gzip.open(file_path, "rt") as f:
        for i, line in enumerate(f):
            # In FASTQ, the sequence is the 2nd line of every 4-line block
            if i % 4 == 1:
                count += len(line.strip())
    
    print(f"Library {library_name}: {count} bases")
    return count

def to_bam_pdf(path_to_pass, bars:list, sz:str):
    ana_path = f"{path_to_pass}/analysis/"
    if not os.path.exists(os.path.dirname(ana_path)):
        os.makedirs(os.path.dirname(ana_path))
    for bar in bars:
        bar_path = f"{ana_path}/{bar}/"
        if not os.path.exists(os.path.dirname(bar_path)):
            os.makedirs(os.path.dirname(bar_path))
        name = f"b{bar}_{sz}"
        genome = "/Users/kalebzuckerman/Desktop/paper_pushing/GCF_000002765.6_GCA_000002765_genomic.fna"
        # genome = "/nfs/jbailey5/baileyweb/bailey_share/resources/MIP_species_resources/pf/Pf_3D7/genomes/unmasked_genome.fa"
        if not os.path.exists(f"{ana_path}/{bar}/{name}.sorted.bam"):
            call = shlex.split(f"touch {ana_path}/{bar}/{name}.sorted.bam")
            subprocess.run(call)
        print("minimap")
        call = f"minimap2 -a -x map-ont {genome} {path_to_pass}/{name}.fastq.gz | samtools view -S -b | samtools sort > {ana_path}/{bar}/{name}.sorted.bam"
        subprocess.call(call, shell=True)
        print("samtools")
        subprocess.call(f"samtools index {ana_path}/{bar}/{name}.sorted.bam {ana_path}/{bar}/{name}.sorted.bam.bai", shell=True)
        print("qualimap")
        subprocess.call(f"qualimap bamqc -bam {ana_path}/{bar}/{name}.sorted.bam -outfile {name}_report.pdf -outdir {ana_path}/{bar}/{name}.sorted_stats", shell=True)
        with open(f"{ana_path}/{bar}/{name}.sorted_stats/genome_results.txt", "r") as f:
            lines = f.readlines()
            to_write = []
            for line in lines:
                reg1 = r'\d*.\d*%'
                reg2 = r'\d*X'
                perc = re.findall(reg1, line)
                cov = re.findall(reg2, line)
                if perc != [] and cov != []:
                    to_write.append((perc[0].rstrip("%"), cov[0].rstrip("X")))
            with open(f"{ana_path}/{bar}/{name}_covPercs.txt", 'w') as f:
                for tup in to_write:
                    f.write(f"{tup[1]},{tup[0]}\n")

def sam_depth_compressor(f_in:str, f_out:str):
    with open(f_in, 'r') as fin, open(f_out, 'w') as fout:
        line = fin.readline()
        if not line:
            return
            
        first = line.strip().split()
        curr_chr = first[0]
        depth = first[2]
        start = int(first[1])
        end = start + 1
        
        for line in fin:
            line = line.strip().split(sep="\t")
            # Logic: If chromosome changes, or there's a gap, or depth changes
            if line[0] != curr_chr or int(line[1]) > end + 1 or line[2] != depth:
                # 1. Write the interval that just finished using the OLD variables
                fout.write(f"{curr_chr}\t{start}\t{end}\t{depth}\n")
                
                # 2. NOW update variables to the new line's data
                curr_chr = line[0]
                start = int(line[1])
                end = start + 1
                depth = line[2]
            else:
                # Same interval, just extend the end coordinate
                end = int(line[1]) + 1
        
        # 3. Write the final interval remaining after the loop finishes
        fout.write(f"{curr_chr}\t{start}\t{end}\t{depth}\n")

def bamcov(path_to_pass, bars, sz):
    for bar in bars:
        name = f"b{bar}_{sz}"
        # subprocess.call(f"bamCoverage --bam {path_to_pass}/analysis/{bar}/{name}.sorted.bam -o {path_to_pass}/analysis/{bar}/{name}_bin10.sorted.bedgraph --outFileFormat bedgraph --binSize 10 --effectiveGenomeSize 23333868", shell=True)
        # subprocess.call(f"bamCoverage --bam {path_to_pass}/analysis/{bar}/{name}.sorted.bam -o {path_to_pass}/analysis/{bar}/{name}_bin500.sorted.bedgraph --outFileFormat bedgraph --binSize 500 --effectiveGenomeSize 23333868", shell=True)
        # dat = pd.read_csv(f"{path_to_pass}/analysis/{bar}/{name}_bin500.sorted.bedgraph", sep='\t', header=None)
        # dat[3] = log10(dat[3] + 1)
        # dat.to_csv(f"{path_to_pass}/analysis/{bar}/{name}_log_bin500.sorted.bedgraph", sep='\t', header=False, index=False)
        # call = shlex.split(f"samtools depth -a -b /Users/kalebzuckerman/Desktop/paper_pushing/3D7_vars_RefSeq.bed {path_to_pass}/analysis/{bar}/{name}.sorted.bam -o {path_to_pass}/analysis/{bar}/{name}_varCov_long.bed")
        # subprocess.call(call)
        sam_depth_compressor(f"{path_to_pass}/analysis/{bar}/{name}_varCov_long.bed", f"{path_to_pass}/analysis/{bar}/{name}_varCov.bed")
        # call = shlex.split(f"samtools depth -a -b /Users/kalebzuckerman/Desktop/paper_pushing/3D7_DR_genes_RefSeq.bed {path_to_pass}/analysis/{bar}/{name}.sorted.bam -o {path_to_pass}/analysis/{bar}/{name}_DRCov_long.bed")
        # subprocess.call(call)
        sam_depth_compressor(f"{path_to_pass}/analysis/{bar}/{name}_DRCov_long.bed", f"{path_to_pass}/analysis/{bar}/{name}_DRCov.bed")

def rerun(path_to_pass, bars, sz):
    for bar in bars:
        name = f"b{bar}_{sz}"
        print(f"samtools depth -a -b /Users/kalebzuckerman/Desktop/paper_pushing/3D7_vars.bed {path_to_pass}/analysis/{bar}/{name}.sorted.bam -o {path_to_pass}/analysis/{bar}/{name}_varCov_long.bed")
        call = shlex.split(f"samtools depth -a -b /Users/kalebzuckerman/Desktop/paper_pushing/3D7_vars_RefSeq.bed {path_to_pass}/analysis/{bar}/{name}.sorted.bam -o {path_to_pass}/analysis/{bar}/{name}_varCov_long.bed")
        subprocess.call(call)
        sam_depth_compressor(f"{path_to_pass}/analysis/{bar}/{name}_varCov_long.bed", f"{path_to_pass}/analysis/{bar}/{name}_varCov.bed")
        print("second")
        call = shlex.split(f"samtools depth -a -b /Users/kalebzuckerman/Desktop/paper_pushing/3D7_DR_genes_RefSeq.bed {path_to_pass}/analysis/{bar}/{name}.sorted.bam -o {path_to_pass}/analysis/{bar}/{name}_DRCov_long.bed")
        subprocess.call(call)
        sam_depth_compressor(f"{path_to_pass}/analysis/{bar}/{name}_DRCov_long.bed", f"{path_to_pass}/analysis/{bar}/{name}_DRCov.bed")


if __name__ == "__main__":
    path_to_data = sys.argv[1] # The folder containing your merge_ files
    libraries = sys.argv[2:]   # List of names (e.g., 'lib1' for merge_lib1.fastq.gz)

    for lib in libraries:
        print(f"Counting {os.path.basename(lib)}")
        # 1. Count bases for this specific library
        num_bases = counter(path_to_data, lib)
        
        # 2. Determine size suffix (KB, MB, GB)
        digits = len(str(num_bases))
        if digits > 9:
            sz = str(round(num_bases/(10**9), 2)) + "GB"
        elif digits > 6:
            sz = str(round(num_bases/(10**6), 2)) + "MB"
        elif digits > 3:
            sz = str(round(num_bases/(10**3), 2)) + "KB"
        else:
            sz = str(num_bases) + "b"

        # # 3. Create a symbolic link or rename for the analysis step
        # # This keeps the rest of your to_bam_pdf function compatible
        # analysis_ready_name = f"b{lib}_{sz}.fastq.gz"
        # os.link(os.path.join(path_to_data, f"{lib}.fastq.gz"), 
        #         os.path.join(path_to_data, analysis_ready_name))

        # # 4. Run Analysis using the new sz suffix
        # print("to bam...")
        # to_bam_pdf(path_to_data, [lib], sz)
        print("bamcov")
        bamcov(path_to_data, [lib], sz)
