# COATswga (Coverage Optimizing and Accurate Toolkit for selective whole genome amplification)

## Overview

COATswga is a command-line tool for identifying coverage-optimizing RNA primer sets for use in selective whole-genome amplification (sWGA).

The pipeline executes in three main steps: 1) finding all potential k-mers using KMC3, 2) filtering out unwanted primers based on parameters both in the program and specified by the users, and 3) forming sets of primers that theoretically tile the genome without forming primer-primer dimers or self-dimers.

Step 3 uses a highly efficient brute-force approach to solve the tiling problem. The algorithm treats each base in the target genome as a sequential list of indices. It records the position of each possible primer and blocks out chunks of indices as "covered" based on the predicted fragment length of the enzyme being used (specified by the user). It leverages `bedtools` for extremely efficient interval computations, cutting down runtime and memory usage drastically while guaranteeing even coverage.

## Installation
First, clone the repository:
```
git clone https://github.com/bailey-lab/coatswga.git
```
Create a virtual environment with the required dependencies with conda:
```
cd coatswga
conda update conda
conda env create -f environment.yml
conda activate coatswga
```
Install `coatswga` with pip:
```
pip install -e .
```
Test the installation with:
```
coatswga -h
```

## Usage
It is recommended to copy the `params.json` file and fill it out with the desired hyperparameters. Command line options can be used but they are not saved and are far less efficient.

If using the JSON format, running is simply
```
coatswga -j /path/to/params.json
```
Without the JSON file, run with:
```
coatswga [options]
```
> [!NOTE]
> The first run with the same foreground/background and k-mer range will always take the longest as KMC3 has to be run. Subsequent runs will finish quickly, allowing edits of the other parameters to fine-tune the primer set.

The full list of hyperparameters and options is:

| JSON entry/long option | Short option | Default value | Description |
| ---------- | ------------ | ------------- | ----------- |
| ---- | -j | None | Path to the file containing parameters in JSON format |
| data_dir | -d | None | Path to the directory to store kmer and CSV files | 
| write | -w | False | Whether or not CSVs of the primers with counts should be written |
| verbose | -v | False | Whether or not extra output is wanted (intermediate sets, timing of parts) | 
| bg_genomes | -bg | None | File paths of the background (off-target) genomes |
| fg_genomes | -fg | None | File paths of the foreground (target) genomes |
| bg_prefix | -bp | None | Prefix for the background files |
| fg_prefixes | -fg | None | List of prefixes for the foreground files |
| cpus | -c | 1 | Number of CPUs to use |
| min_primer_length | -m | 8 | Minimum primer length |
| max_primer_length | -M | 16 | Maximum primer length |
| min_tm | -t | 15 | Minimum predicted primer melting temperature |
| max_tm | -T | 45 | Maximum predicted primer melting temperature | 
| min_fg_count | -g | 200 | Minimum number of foreground genome occurrences, helps to filter out less effective primers |
| max_ratio | -r | 1 | Maximum value of the background occurrences divided by foreground occurrences. Lower values mean higher specificity of the primers but also cuts out potentially valuable less-specific primers. | 
| fragment_length | -l | 10000 | Theoretical length of DNA fragments formed by DNA polymerase enzymes |
| target_coverage | -o | 0.95 | Proportion of the genome to cover |
| force_coverage_threshold | -f | False | Force the program to run until either every potential primer has been checked or the target coverage threshold has been reached. Leaving this as false allows the program to exit after only the first several primers have been checked even if the target coverage threshold has not been reached |
| min_set_size | -s | 10 | Minimum number of primers to include in the set |
| existing_primers | -p | [ ] | Any primers currently in use that the program should check for potential primer-primer dimers |
 
The JSON file can be used in tandem with command line options as long as `-j` is specified. Command line options will override but not overwrite anything in the JSON file.

> [!NOTE]
> `min_fg_count`, `max_ratio`, and `force_coverage_threshold` all have large impacts on the runtime in exchange for specificity/coverage.
>
> `min_fg_count` puts a floor on how many foreground occurrences are needed for the primer to pass the filters. Raising this value reduces the number of primers to check thereby reducing runtime, but also potentially eliminates primers needed to fill gaps in the coverage of others.
>
> `max_ratio` puts a ceiling on how many background hits per foreground hit a primer can have. Lowering this value guarantees increased specificity and a lower runtime but risks losing primers needed to cover larger swaths of the target genome
>
> `force_coverage_threshold` is only recommended as a last-resort option as it can make the program run for days. If multiple runs have been made while tuning the other hyperparameters without satisfying results, this option can be made `True` to force the program to continue to evaluate primer sets until one is found that covers the target threshold. When `False`, the program only evaluates the number of top primers equal to the number of CPUs, for the sake of efficiency.
>
> `existing_primers` is meant to be used as a hole-filling option. If there is already a primer set in use that does poorly on a couple smaller areas, those sequences can be passed in with the primer set specified at this option and the program with ensure any new primers given will not form dimers with any primers in the set

## Resources
COATswga uses the following modules:
* `KMC3`. K-mer counting and filtering. https://github.com/refresh-bio/KMC
* `bedtools`. Interval arithmetic. https://bedtools.readthedocs.io/en/latest/

Primer-dimer calculations were based on an algorithm implemented by Jason Hendry in `multiply` https://github.com/JasonAHendry/multiply/tree/master
