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

def counter(path_to_barcode:str, bar, num_bases:int, sz:str):
    # print(f"Bar: {bar}")
    if len(bar) == 2:
        path = path_to_barcode + "barcode" + bar
    else:
        path = path_to_barcode + bar
    # print(path)
    if num_bases != 0:
        if os.path.exists(f"{path}/b{bar}_{sz}.fastq.gz"):
            call = shlex.split(f"rm {path}/b{bar}_{sz}.fastq.gz")
            subprocess.call(call)
        call = shlex.split(f"touch {path}/b{bar}_{sz}.fastq.gz")
        subprocess.run(call)
    id_search = r'fastq_runid_([a-zA-Z0-9]+)_(\d+)_\d*.fastq.gz'
    fn = 0
    count = 0
    zips = []
    files = os.listdir(path)
    for file in files:
        # print(file)
        result = re.search(id_search, file)
        if result:
            # id = result.group(1)
            # fn = result.group(2)
            # print(id)
            # file = f"fastq_runid_{id}_{fn}_0.fastq.gz"
            # call = shlex.split(f"gunzip {path}/{file}")
            # subprocess.call(call)
            with gzip.open(f"{path}/{file}", "rt") as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    line = line.strip()
                    if "read" in line:
                        count += len(lines[i + 1])
            # call = shlex.split(f"gzip {path}/{file.rstrip('.gz')}")
            # subprocess.call(call)
            zips.append(f"{path}/{file}")
            # print(f"\rfile: {file}   count: {count}", end=" ")
            if num_bases != 0:
                subprocess.run(f"cat {path}/{file} >> {path}/b{bar}_{sz}.fastq.gz", shell=True, stdout=subprocess.DEVNULL)
            if num_bases > 0 and count >= num_bases:
                print("")
                break
    print(f"{count} bases")
    return count

def to_bam_pdf(path_to_pass, bar, sz:str):
    ana_path = f"{path_to_pass}/stepped_analysis"
    if not os.path.exists(os.path.dirname(ana_path)):
        os.makedirs(os.path.dirname(ana_path))
    if len(bar) == 2:
        named_bar = f"barcode{bar}"
    else:
        named_bar = f"{bar}"
    bar_path = f"{ana_path}/{named_bar}/"
    if not os.path.exists(os.path.dirname(bar_path)):
        os.makedirs(os.path.dirname(bar_path))
    name = f"b{bar}_{sz}"
    # genome = "/Users/kaleb/Desktop/Bailey_Lab/code/newswga/genomes/malaria_GCF_000002765.5_GCA_000002765_genomic.fna"
    genome = "/nfs/jbailey5/baileyweb/bailey_share/resources/MIP_species_resources/pf/Pf_3D7/genomes/unmasked_genome.fa"
    if not os.path.exists(f"{bar_path}/{name}.sorted.bam"):
        call = shlex.split(f"touch {bar_path}/{name}.sorted.bam")
        subprocess.run(call)
    print("minimap")
    call = f"minimap2 -a -x map-ont {genome} {path_to_pass}/{named_bar}/{name}.fastq.gz | samtools view -S -b | samtools sort > {bar_path}/{name}.sorted.bam"
    subprocess.call(call, shell=True)
    print("samtools")
    subprocess.call(f"samtools index {bar_path}/{name}.sorted.bam {bar_path}/{name}.sorted.bam.bai", shell=True)
    print("qualimap")
    subprocess.call(f"qualimap bamqc -bam {bar_path}/{name}.sorted.bam -outfile {name}_report.pdf -outdir {bar_path}/{name}.sorted_stats", shell=True)
    with open(f"{bar_path}{name}.sorted_stats/genome_results.txt", "r") as f:
            lines = f.readlines()
            to_write = []
            for line in lines:
                reg1 = r'\d*.\d*%'
                reg2 = r'\d*X'
                perc = re.findall(reg1, line)
                cov = re.findall(reg2, line)
                if perc != [] and cov != []:
                    to_write.append((perc[0].rstrip("%"), cov[0].rstrip("X")))
            with open(f"{bar_path}{name}_covPercs.txt", 'w') as f:
                for tup in to_write:
                    f.write(f"{tup[1]},{tup[0]}\n")

if __name__ == "__main__":
    paths = [
        "/nfs/jbailey5/baileyweb/kzucker2/sWGA_seq/combine_0831_1116_analysis/", 
        "/nfs/jbailey5/baileyweb/kzucker2/sWGA_seq/231117_sWGA_P2/pass/",
        "/nfs/jbailey5/baileyweb/kzucker2/sWGA_seq/230831_sWGA/pass/",

    ]
    all_bars = [
        ["coat", "soap", "sanger"],
        ["17", "19", "24"],
        ["19", "22", "23"]
    ]

    bases = 350000000
    
    coverages = [1, 10]

    for i in range(len(paths)):
        path_to_pass = paths[i]
        bars = all_bars[i]
        if path_to_pass[-1] == "/":
            path = path_to_pass
        else:
            path = path_to_pass + "/"
        total_bases = []
        for bar in bars:
            total_bases.append(counter(path, bar, 0, ""))
        bar_to_covs = {bar:[] for bar in bars}
        for i, count in enumerate(total_bases):
            bar = bars[i]
            covs = {cover:[] for cover in coverages}
            print(f"Bar: {bar}")
            if bases <= count:
                digits = len(str(bases))
                if digits > 9:
                    sz = str(int(bases/(10**9))) + "GB"
                elif digits > 6:
                    sz = str(int(bases/(10**6))) + "MB"
                elif digits > 3:
                    sz = str(int(bases/(10**3))) + "KB"
                else:
                    sz = str(bases) + "B"
                print(sz)
                counter(path, bar, bases, sz)
                to_bam_pdf(path_to_pass, bar, sz)
                if len(bar) == 2:
                    bar_path = f"{path}stepped_analysis/barcode{bar}/"
                else:
                    bar_path = f"{path}stepped_analysis/{bar}/"
                df = pd.read_csv(f"{bar_path}b{bar}_{sz}_covPercs.txt", sep=",", header=None, index_col=None)
                for coverage in coverages:
                    covs[coverage].append((int(bases/1000000),df.iloc[coverage - 1, 1]))
            for coverage in coverages:
                with open(f"{path}b{bar}_{coverage}X_stepped_cov.txt", 'w') as f:
                    for cov in covs[coverage]:
                        f.write(f"{cov[0]},{cov[1]}\n")