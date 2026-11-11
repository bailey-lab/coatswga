import pandas as pd

path = "/Users/kalebzuckerman/Desktop/paper_pushing/combined_analysis/stepped"
bars = ["coat", "soap", "sanger"]

counts = [10, 20, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]
labs = [str(a) + "MB" for a in [10, 20, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900]] + ["1GB"]
coverages = [1, 5, 10]
for bar in bars:
    for coverage in coverages:
        to_write = ""
        for count, lab in zip(counts, labs):
            df = pd.read_csv(f"{path}/b{bar}_{lab}_covPercs.txt", sep=",", header=None, index_col=None)
            to_write += f"{count},{df.iloc[coverage - 1, 1]}\n"
        with open(f"{path}/{bar}_{coverage}X_stepped_cov.csv", "w") as f:
            f.write(to_write)