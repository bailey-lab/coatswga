import os
import sys
import re
import subprocess
import shlex
from numpy import log10
import gzip

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
    # fin = sys.argv[1]
    # fout = sys.argv[2]
    # sam_depth_compressor(fin, fout)
    dir = "/Users/kalebzuckerman/Desktop/paper_pushing/combined_analysis/longs"
    for f_name in os.listdir(dir):
        if f_name[0] != ".":
            out = f_name.replace("_long", "")
            sam_depth_compressor(os.path.join(dir, f_name), os.path.join(dir, out))