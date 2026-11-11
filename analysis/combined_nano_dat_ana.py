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

def counter(path_to_barcode:str, bars:list, num_bases:int, sz:str):
    total_bases = []
    for bar in bars:
        print(f"Bar: {bar}")
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
        total_bases.append(count)
    return total_bases

def to_bam_pdf(path_to_pass, bars:list, sz:str):
    ana_path = f"{path_to_pass}/analysis/"
    if not os.path.exists(os.path.dirname(ana_path)):
        os.makedirs(os.path.dirname(ana_path))
    for bar in bars:
        bar_path = f"{ana_path}/{bar}/"
        if not os.path.exists(os.path.dirname(bar_path)):
            os.makedirs(os.path.dirname(bar_path))
        name = f"b{bar}_{sz}"
        # genome = "/Users/kaleb/Desktop/Bailey_Lab/code/newswga/genomes/malaria_GCF_000002765.5_GCA_000002765_genomic.fna"
        genome = "/nfs/jbailey5/baileyweb/bailey_share/resources/MIP_species_resources/pf/Pf_3D7/genomes/unmasked_genome.fa"
        if not os.path.exists(f"{ana_path}/{bar}/{name}.sorted.bam"):
            call = shlex.split(f"touch {ana_path}/{bar}/{name}.sorted.bam")
            subprocess.run(call)
        print("minimap")
        call = f"minimap2 -a -x map-ont {genome} {path_to_pass}/{bar}/{name}.fastq.gz | samtools view -S -b | samtools sort > {ana_path}/{bar}/{name}.sorted.bam"
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
        first = fin.readline().strip().split()
        chr = first[0]
        depth = first[2]
        start = int(first[1])
        end = start + 1
        for line in fin:
            line = line.strip().split(sep="\t")
            if line[0] != chr or int(line[1]) > end + 1 or line[2] != depth:
                fout.write(f"{line[0]}\t{start}\t{end}\t{depth}\n")
                start = int(line[1])
                end = start + 1
                chr = line[0]
                depth = line[2]
            else:
                end = int(line[1]) + 1
        fout.write(f"{line[0]}\t{start}\t{end}\t{depth}\n")

def bamcov(path_to_pass, bars, sz):
    for bar in bars:
        name = f"b{bar}_{sz}"
        subprocess.call(f"bamCoverage --bam {path_to_pass}/analysis/{bar}/{name}.sorted.bam -o {path_to_pass}/analysis/{bar}/{name}_bin10.sorted.bedgraph --outFileFormat bedgraph --binSize 10 --effectiveGenomeSize 23333868", shell=True)
        subprocess.call(f"bamCoverage --bam {path_to_pass}/analysis/{bar}/{name}.sorted.bam -o {path_to_pass}/analysis/{bar}/{name}_bin500.sorted.bedgraph --outFileFormat bedgraph --binSize 500 --effectiveGenomeSize 23333868", shell=True)
        dat = pd.read_csv(f"{path_to_pass}/analysis/{bar}/{name}_bin500.sorted.bedgraph", sep='\t', header=None)
        dat[3] = log10(dat[3])
        dat.to_csv(f"{path_to_pass}/analysis/{bar}/{name}_log_bin500.sorted.bedgraph", sep='\t', header=False, index=False)
        call = shlex.split(f"samtools depth -a -b /nfs/jbailey5/baileyweb/kzucker2/3D7_vars.bed {path_to_pass}/analysis/{bar}/{name}.sorted.bam -o {path_to_pass}/analysis/{bar}/{name}_varCov_long.bed")
        subprocess.call(call)
        sam_depth_compressor(f"{path_to_pass}/analysis/{bar}/{name}_varCov_long.bed", f"{path_to_pass}/analysis/{bar}/{name}_varCov.bed")
        call = shlex.split(f"samtools depth -a -b /nfs/jbailey5/baileyweb/kzucker2/3D7_DR_genes.bed {path_to_pass}/analysis/{bar}/{name}.sorted.bam -o {path_to_pass}/analysis/{bar}/{name}_DRCov_long.bed")
        subprocess.call(call)
        sam_depth_compressor(f"{path_to_pass}/analysis/{bar}/{name}_DRCov_long.bed", f"{path_to_pass}/analysis/{bar}/{name}_DRCov.bed")

def clair3(path_to_pass, bars, sz):
    out_path = f"{path_to_pass}/clair_out"
    if not os.path.exists(os.path.dirname(out_path)):
        os.makedirs(os.path.dirname(out_path))
    
    pass

if __name__ == "__main__":
    path_to_pass = sys.argv[1]
    path = path_to_pass + "/"
    bars = [bar for bar in sys.argv[2:]]
    print(bars)
    total_bases = counter(path, bars, 0, "")
    least = min(total_bases)
    digits = len(str(least))
    if digits > 9:
        sz = str(int(least/(10**9))) + "GB"
    elif digits > 6:
        sz = str(int(least/(10**6))) + "MB"
    elif digits > 3:
        sz = str(int(least/(10**3))) + "KB"
    else:
        sz = str(least) + "b"
    print(sz)
    counter(path, bars, least, sz)
    to_bam_pdf(path_to_pass, bars, sz)
    bamcov(path_to_pass, bars, sz)
    
# 0%   10   20   30   40   50   60   70   80   90   100%
# |----|----|----|----|----|----|----|----|----|----|