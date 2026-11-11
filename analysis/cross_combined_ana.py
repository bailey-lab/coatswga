import os
import sys
import re
import subprocess
import shlex
import pandas as pd
from numpy import log10
import gzip

# Needs module load minimap2/2.24 samtools/1.9 qualimap/2.2.1 java/jdk-11 R/3.5.2
# also needs conda environment with bamcov (deeptools) and clair3

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
                fout.write(f"{chr}\t{start}\t{end}\t{depth}\n")
                start = int(line[1])
                end = start + 1
                chr = line[0]
                depth = line[2]
            else:
                end = int(line[1]) + 1
        fout.write(f"{line[0]}\t{start}\t{end}\t{depth}\n")

if __name__ == "__main__":
    comb_dir = "/nfs/jbailey5/baileyweb/kzucker2/sWGA_seq/combine_0831_1116_analysis/cross_combine_analysis/"
    bams = ["/nfs/jbailey5/baileyweb/kzucker2/sWGA_seq/combine_0831_1116_analysis/analysis/coat/bcoat_2GB.sorted.bam",
            comb_dir + "coatxlong_2GB.sorted.bam", comb_dir + "coatxsanger_2GB.sorted.bam", comb_dir + "coatxsoap_2GB.sorted.bam"]
    ana_path = comb_dir + "analysis"
    for bam, name in zip(bams, ["coat_2GB", "coatxlong_2GB", "coatxsanger_2GB", "coatxsoap_2GB"]):
        subprocess.call(f"qualimap bamqc -bam {bam} -outfile {name}_report.pdf -outdir {ana_path}/{name}.sorted_stats", shell=True)
        with open(f"{ana_path}/{name}.sorted_stats/genome_results.txt", "r") as f:
            lines = f.readlines()
            to_write = []
            for line in lines:
                reg1 = r'\d*.\d*%'
                reg2 = r'\d*X'
                perc = re.findall(reg1, line)
                cov = re.findall(reg2, line)
                if perc != [] and cov != []:
                    to_write.append((perc[0].rstrip("%"), cov[0].rstrip("X")))
            with open(f"{ana_path}/{name}_covPercs.txt", 'w') as f:
                for tup in to_write:
                    f.write(f"{tup[1]},{tup[0]}\n")
        subprocess.call(f"bamCoverage --bam {bam} -o {ana_path}/{name}_bin10.sorted.bedgraph --outFileFormat bedgraph --binSize 10 --effectiveGenomeSize 23333868", shell=True)
        call = shlex.split(f"samtools depth -a -b /nfs/jbailey5/baileyweb/kzucker2/3D7_vars.bed {bam} -o {ana_path}/{name}_varCov_long.bed")
        subprocess.call(call)
        sam_depth_compressor(f"{ana_path}/{name}_varCov_long.bed", f"{ana_path}/{name}_varCov.bed")
        call = shlex.split(f"samtools depth -a -b /nfs/jbailey5/baileyweb/kzucker2/3D7_DR_genes.bed {bam} -o {ana_path}/{name}_DRCov_long.bed")
        subprocess.call(call)
        sam_depth_compressor(f"{ana_path}/{name}_DRCov_long.bed", f"{ana_path}/{name}_DRCov.bed")

    
# 0%   10   20   30   40   50   60   70   80   90   100%
# |----|----|----|----|----|----|----|----|----|----|