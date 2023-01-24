# ggCaller: a bacterial gene caller for pangenome graphs <img src='docs/images/ggCaller_logo.png' align="right" height="100" />

ggCaller traverses [Bifrost](https://github.com/pmelsted/bifrost) graphs constructed from bacterial genomes to identify putative gene sequences, known as open reading frames (ORFs). 

ggCaller incorporates [Balrog](https://github.com/salzberg-lab/Balrog) to filter ORFs to improve specificity of calls and [Panaroo](https://github.com/gtonkinhill/panaroo) for pangenome analysis and quality control.

## Documentation

Guides for installation, usage and a tutorial can be found [here](https://ggcaller.readthedocs.io/en/latest/).

## Installation

ggCaller is available on Linux. If you are running Windows 10/11, Linux can be installed via the Windows Subsystem for Linux ([WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10)).

We plan to get a MacOS version up and running in the future.

NOTE: We are aware of issues installing from conda at the moment. We recommend installing from source at this time.

### Installation from source (recommended)
Required packages and versions can be found in ```environment_linux.yml``` and ```environment_macOS.yml``` depending on your operating system. In addition, a C++17 compiler (e.g. gcc >=7.3) is required.

For example, using conda (creates ```ggc_env``` environment)

```
conda env create -f environment_linux.yml
conda activate ggc_env
```

Once all required packages are installed, install ggCaller using:

```
git clone --recursive https://github.com/samhorsfield96/ggCaller
cd ggCaller
python setup.py install
```

### Installation via conda

Install through [bioconda](http://bioconda.github.io/):

```conda install ggcaller```

If conda is not installed, first install [miniconda](https://docs.conda.io/en/latest/miniconda.html), then add the correct channels:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

## Citation

If you use this code, please cite:

Bifrost: 
Holley, G., Melsted, P. Bifrost: highly parallel construction and indexing of colored and compacted de Bruijn graphs. Genome Biol 21, 249 (2020). https://doi.org/10.1186/s13059-020-02135-8

Balrog:
Sommer MJ, Salzberg SL. Balrog: A universal protein model for prokaryotic gene prediction. PLoS Comput Biol 17(2): e1008727 (2021). https://doi.org/10.1371/journal.pcbi.1008727

Panaroo:
Tonkin-Hill, G., MacAlasdair, N., Ruis, C. et al. Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biol 21, 180 (2020). https://doi.org/10.1186/s13059-020-02090-4

SDSL v3:
[Succinct Data Structure Library 3.0](https://github.com/xxsds/sdsl-lite)

Edlib:
Šošić, M., Šikić, M. Edlib: a C/C++ library for fast, exact sequence alignment using edit distance. Bioinformatics 33, 9, (2017). https://doi.org/10.1093/bioinformatics/btw753

Eigen v3:
Guennebaud, G., Jacob, B. et al. Eigen v3 (2010). http://eigen.tuxfamily.org