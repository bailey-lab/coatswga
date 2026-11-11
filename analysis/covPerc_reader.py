import pandas as pd
import sys


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

steps = [1000000 * x for x in [10, 20, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]]

coverages = [1, 10]
for i in range(len(paths)):
    path = paths[i]
    bars = all_bars[i]
    for bar in bars:
        if len(bar) == 2:
            bar_seg = f"barcode{bar}"
        else:
            bar_seg = bar
        covs = {cover:[] for cover in coverages}
        for bases in steps:
            digits = len(str(bases))
            if digits > 9:
                sz = str(int(bases/(10**9))) + "GB"
            elif digits > 6:
                sz = str(int(bases/(10**6))) + "MB"
            elif digits > 3:
                sz = str(int(bases/(10**3))) + "KB"
            else:
                sz = str(bases) + "B"
            try:
                df = pd.read_csv(f"{path}stepped_analysis/{bar_seg}/b{bar}_{sz}_covPercs.txt", sep=",", header=None, index_col=None)
                for coverage in coverages:
                    covs[coverage].append((int(bases/1000000),df.iloc[coverage - 1, 1]))
            except:
                pass
        for coverage in coverages:
            with open(f"{path}b{bar}_{coverage}X_stepped_cov.txt", 'w') as f:
                for cov in covs[coverage]:
                    f.write(f"{cov[0]},{cov[1]}\n")
