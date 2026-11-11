import os

path = "/nfs/jbailey5/baileyweb/kzucker2/sWGA_seq/230831_sWGA/pass/barcode23"
# path = "/Users/kaleb/Desktop/Bailey_Lab/code"
id = "7f409fa7d63c15f0de3bf5913bc1a8d5e10b7de3"

n_bases = 360000000

fn = 0
count = 0
zips = []
while os.path.exists(f"{path}/fastq_runid_{id}_{fn}_0.fastq.gz"):
    file = f"fastq_runid_{id}_{fn}_0.fastq.gz"
    os.system(f"gunzip {path}/{file}")
    with open(f"{path}/{file.rstrip('.gz')}", "r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            line = line.strip()
            if "read" in line:
                count += len(lines[i + 1])
    os.system(f"gzip {path}/{file.rstrip('.gz')}")
    zips.append(file)
    print(f"\rfile: {file}   count: {count}", end=" ")
    if count >= n_bases:
        print("")
        break
    fn += 1
tot = ""
for z in zips:
    tot = tot + f"{path}/{z} "
print("pre cat")
print(tot)
os.system(f"cat {tot.strip()} > {path}/bases_{n_bases}.fastq.gz")
