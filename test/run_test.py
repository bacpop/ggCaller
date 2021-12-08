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
sys.stderr.write("Running reference build workflow\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --max-path-length 5000 --clean-mode strict --min-orf-length 100 --max-ORF-overlap 55 --alignment core --aligner def --annotation fast --evalue 0.0001 --search_radius 3000 --save",
    shell=True,
    check=True)

sys.stderr.write("Running unitig query workflow\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --data test_dir/ggc_data --query pneumo_CL_group2_queries.fasta --out test_dir",
    shell=True,
    check=True)

### reference read workflow ###
sys.stderr.write("Running reference read workflow\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --len-diff-cutoff 0.99 --alignment pan --aligner def --annotation sensitive --evalue 0.01 --high_var_flag 4 --core_threshold 0.96",
    shell=True,
    check=True)

### reference read workflow ###
sys.stderr.write("Running reference read workflow with custom annotation dbs\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --diamonddb ../db/diamond/Uniprot_test.fasta --hmmdb ../db/hmm/HAMAP.hmm --annotation sensitive",
    shell=True,
    check=True)

### read build workflow ###
sys.stderr.write("Running read build workflow\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --min-orf-score 150 --max-orf-orf-distance 5000 --identity-cutoff 0.99 --clean-mode moderate --length_outlier_support_proportion 0.15 --min_edge_support_sv 1",
    shell=True,
    check=True)

### reads + refs build workflow ###
sys.stderr.write("Running reads + reference build workflow\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --out test_dir --repeat --family_threshold 0.75 --merge_paralogs --refind_prop_match 0.15 --edge_support_threshold 2 --remove_by_consensus True",
    shell=True,
    check=True)

### reads read workflow ###
sys.stderr.write("Running reads read workflow\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir --repeat --clean-mode sensitive --min_trailing_support 1 --trailing_recursive 1",
    shell=True,
    check=True)

sys.stderr.write("Running reference build without panaroo\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --max-path-length 5000 --min-orf-length 100 --repeat --no-clustering",
    shell=True,
    check=True)

### reference read workflow ###
sys.stderr.write("Running reference read workflow without balrog\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --no-filter --len-diff-cutoff 0.99 --alignment pan --aligner def --annotation sensitive --truncation_threshold 0.7 --ignore_pseduogenes --no_variants",
    shell=True,
    check=True)

sys.stderr.write("Tests completed\n")
