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

### reference build workflows ###
sys.stderr.write("Running reference build workflow without alignment\n")
subprocess.run(python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir", shell=True,
               check=True)

sys.stderr.write("Running reference build workflow with default pangenome alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner def",
    shell=True, check=True)

sys.stderr.write("Running reference build workflow with reference-guided pangenome alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner ref",
    shell=True, check=True)

sys.stderr.write("Running reference build workflow with default core alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner def",
    shell=True, check=True)

sys.stderr.write("Running reference build workflow with reference-guided core alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner ref",
    shell=True, check=True)

# fast annotation
sys.stderr.write("Running reference build workflow with default pangenome alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner def --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running reference build workflow with reference-guided pangenome and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner ref --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running reference build workflow with default core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner def --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running reference build workflow with reference-guided core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner ref --annotation fast",
    shell=True, check=True)

# sensitive annotation
sys.stderr.write("Running reference build workflow with default pangenome alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner def --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running reference build workflow with reference-guided pangenome and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner ref --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running reference build workflow with default core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner def --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running reference build workflow with reference-guided core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner ref --annotation sensitive",
    shell=True, check=True)

### reference read workflows ###
sys.stderr.write("Running reference read workflow without alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir",
    shell=True, check=True)

sys.stderr.write("Running reference read workflow with default pangenome alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --alignment pan --aligner def",
    shell=True, check=True)

sys.stderr.write("Running reference read workflow with reference-guided pangenome alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --alignment pan --aligner ref",
    shell=True, check=True)

sys.stderr.write("Running reference read workflow with default core alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --alignment core --aligner def",
    shell=True, check=True)

sys.stderr.write("Running reference read workflow with reference-guided core alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --alignment core --aligner ref",
    shell=True, check=True)

# fast annotation
sys.stderr.write("Running reference read workflow with default pangenome alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --alignment pan --aligner def --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running reference read workflow with reference-guided pangenome and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --alignment pan --aligner ref --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running reference read workflow with default core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --alignment core --aligner def --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running reference read workflow with reference-guided core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --alignment core --aligner ref --annotation fast",
    shell=True, check=True)

# sensitive annotation
sys.stderr.write("Running reference read workflow with default pangenome alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --alignment pan --aligner def --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running reference read workflow with reference-guided pangenome and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --alignment pan --aligner ref --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running reference read workflow with default core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --alignment core --aligner def --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running reference read workflow with reference-guided core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --out test_dir --alignment core --aligner ref --annotation sensitive",
    shell=True, check=True)

### reads build workflows ###
sys.stderr.write("Running reads build workflow without alignment\n")
subprocess.run(python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir", shell=True,
               check=True)

sys.stderr.write("Running reads build workflow with default pangenome alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner def",
    shell=True, check=True)

sys.stderr.write("Running reads build workflow with reference-guided pangenome alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner ref",
    shell=True, check=True)

sys.stderr.write("Running reads build workflow with default core alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner def",
    shell=True, check=True)

sys.stderr.write("Running reads build workflow with reference-guided core alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner ref",
    shell=True, check=True)

# fast annotation
sys.stderr.write("Running reads build workflow with default pangenome alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner def --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running reads build workflow with reference-guided pangenome and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner ref --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running reads build workflow with default core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner def --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running reads build workflow with reference-guided core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner ref --annotation fast",
    shell=True, check=True)

# sensitive annotation
sys.stderr.write("Running reads build workflow with default pangenome alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner def --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running reads build workflow with reference-guided pangenome and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner ref --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running reads build workflow with default core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner def --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running reads build workflow with reference-guided core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner ref --annotation sensitive",
    shell=True, check=True)

### reads read workflows ###
sys.stderr.write("Running non-reference read workflow without alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir",
    shell=True, check=True)

sys.stderr.write("Running non-reference read workflow with default pangenome alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir --alignment pan --aligner def",
    shell=True, check=True)

sys.stderr.write("Running non-reference read workflow with reference-guided pangenome alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir --alignment pan --aligner ref",
    shell=True, check=True)

sys.stderr.write("Running non-reference read workflow with default core alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir --alignment core --aligner def",
    shell=True, check=True)

sys.stderr.write("Running non-reference read workflow with reference-guided core alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir --alignment core --aligner ref",
    shell=True, check=True)

# fast annotation
sys.stderr.write("Running non-reference read workflow with default pangenome alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir --alignment pan --aligner def --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running non-reference read workflow with reference-guided pangenome and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir --alignment pan --aligner ref --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running non-reference read workflow with default core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir --alignment core --aligner def --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running non-reference read workflow with reference-guided core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir --alignment core --aligner ref --annotation fast",
    shell=True, check=True)

# sensitive annotation
sys.stderr.write("Running non-reference read workflow with default pangenome alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir --alignment pan --aligner def --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running non-reference read workflow with reference-guided pangenome and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir --alignment pan --aligner ref --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running non-reference read workflow with default core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir --alignment core --aligner def --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running non-reference read workflow with reference-guided core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --graph pneumo_CL_group2.gfa --colours pneumo_CL_group2.bfg_colors --not-ref --out test_dir --alignment core --aligner ref --annotation sensitive",
    shell=True, check=True)

### reference + read build workflows ###
sys.stderr.write("Running reference + read build workflow without alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --kmer 31 --out test_dir",
    shell=True, check=True)

sys.stderr.write("Running reference + read build workflow with default pangenome alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner def",
    shell=True, check=True)

sys.stderr.write("Running reference + read build workflow with reference-guided pangenome alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner ref",
    shell=True, check=True)

sys.stderr.write("Running reference + read build workflow with default core alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner def",
    shell=True, check=True)

sys.stderr.write("Running reference + read build workflow with reference-guided core alignment\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner ref",
    shell=True, check=True)

# fast annotation
sys.stderr.write("Running reference + read build workflow with default pangenome alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner def --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running reference + read build workflow with reference-guided pangenome and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner ref --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running reference + read build workflow with default core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner def --annotation fast",
    shell=True, check=True)

sys.stderr.write("Running reference + read build workflow with reference-guided core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner ref --annotation fast",
    shell=True, check=True)

# sensitive annotation
sys.stderr.write("Running reference + read build workflow with default pangenome alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner def --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running reference + read build workflow with reference-guided pangenome and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment pan --aligner ref --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running reference + read build workflow with default core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner def --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Running reference + read build workflow with reference-guided core alignment and fast annotation\n")
subprocess.run(
    python_cmd + " ../ggcaller-runner.py --refs pneumo_CL_group2.txt --reads pneumo_CL_group2.txt --kmer 31 --out test_dir --alignment core --aligner ref --annotation sensitive",
    shell=True, check=True)

sys.stderr.write("Tests completed\n")
