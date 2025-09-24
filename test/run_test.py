#!/usr/bin/env python
# Copyright 2021 Samuel Horsfield

"""Tests for ggCaller"""

import subprocess
import os
import sys
import shutil

if os.environ.get("GGCALLER_PYTHON"):
    python_cmd = os.environ.get("GGCALLER_PYTHON")
else:
    python_cmd = "python"

### reference build workflow ###
sys.stderr.write("Running reference build workflow without annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --max-path-length 5000 --clean-mode strict --min-orf-length 100 --max-ORF-overlap 55 --alignment core --aligner def --annotation none --evalue 0.0001 --search-radius 3000",
    shell=True,
    check=True)

sys.stderr.write("Running reference build workflow generating gff only\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 13 --out test_dir --max-path-length 5000 --clean-mode strict --min-orf-length 0 --max-ORF-overlap 55 --alignment core --aligner def --annotation fast --evalue 0.0001 --search-radius 3000 --save --gene-finding-only",
    shell=True,
    check=True)

sys.stderr.write("Running reference build workflow with annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 13 --out test_dir --max-path-length 5000 --clean-mode strict --min-orf-length 0 --max-ORF-overlap 55 --alignment core --aligner def --annotation fast --evalue 0.0001 --search-radius 3000 --save",
    shell=True,
    check=True)

print("CWD:", os.getcwd())
print("Contents:", os.listdir())
print("test_dir contents:", os.listdir("test_dir"))

# Determine which color file was generated
color_file_candidates = ["pneumo_CL_group2.color.bfg", "pneumo_CL_group2.bfg_colors"]
for cf in color_file_candidates:
    if os.path.isfile(cf):
        color_file = cf
        break
else:
    raise FileNotFoundError(
        "Neither .color.bfg nor .bfg_colors file was found in test_dir"
    )

sys.stderr.write("Running unitig query workflow\n")
subprocess.run(
    python_cmd + f" ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours {color_file} --prev-run test_dir --query pneumo_CL_group2_queries.fasta --query-id 0.5 --out test_dir",
    shell=True,
    check=True)

### reference read workflow ###
sys.stderr.write("Running reference read workflow\n")
subprocess.run(
    python_cmd + f" ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours {color_file} --out test_dir --len-diff-cutoff 0.99 --alignment pan --aligner def --annotation sensitive --evalue 0.01 --core-threshold 0.96",
    shell=True,
    check=True)

### read build workflow ###
sys.stderr.write("Running read build workflow\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --min-orf-score 150 --max-orf-orf-distance 5000 --identity-cutoff 0.99 --annotation ultrasensitive --clean-mode moderate --length-outlier-support-proportion 0.15 --min-edge-support-sv 1",
    shell=True,
    check=True)

### reads + refs build workflow ###
sys.stderr.write("Running reads + reference build workflow\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --out test_dir --repeat --family-threshold 0.75 --merge-paralogs --refind-prop-match 0.15 --edge-support-threshold 2",
    shell=True,
    check=True)

sys.stderr.write("Running reads + reference build workflow generating gff only\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --out test_dir --repeat --family-threshold 0.75 --merge-paralogs --refind-prop-match 0.15 --gene-finding-only",
    shell=True,
    check=True)

### reads read workflow ###
sys.stderr.write("Running reads read workflow\n")
subprocess.run(
    python_cmd + f" ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours {color_file} --not-ref --out test_dir --repeat --clean-mode sensitive --min-trailing-support 1 --trailing-recursive 1",
    shell=True,
    check=True)

sys.stderr.write("Running reference build without panaroo\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --max-path-length 5000 --min-orf-length 100 --repeat --no-clustering",
    shell=True,
    check=True)

### reference read workflow ###
sys.stderr.write("Running reference read workflow without models\n")
subprocess.run(
    python_cmd + f" ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours {color_file} --out test_dir --no-filter --len-diff-cutoff 0.99 --alignment pan --aligner def --annotation sensitive --truncation-threshold 0.7 --ignore-pseduogenes --no-variants",
    shell=True,
    check=True)

sys.stderr.write("Tests completed\n")
