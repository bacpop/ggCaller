# ggCaller: a graph gene caller for Bifrost graphs

Traverses Bifrost graphs to identify putative protein coding sequences, known as open reading frames (ORFs).

ggCaller uses pyGFA to convert Bifrost graphs to Networkx graph objects to enable graph traversal.

## Requirements
- python==3.8
- biopython
- networkx==1.11

## Installation

###Cloning the repository

```git clone --recursive https://github.com/samhorsfield96/ggCaller```

## Input files for ggCaller

ggCaller takes a Bifrost GFA file generated by ```Bifrost build```, and a TSV file generated by ```Bifrost query```. See the Bifrost repository for installation (https://github.com/pmelsted/bifrost).

### Generate Bifrost GFA

Using ```Bifrost build```, generate a GFA file with an associated colours file using a list of reference genomes. Here k-mer size is set to 31 bp by ```-k 31```.

```Bifrost build -k 31 -c -r references.txt -o graph```

### Generate TSV file

The TSV file should be generated by querying the constituent unitigs of the GFA file against the GFA itself.
The accompanying file 'gfa_to_fasta.py' can be used to convert a GFA file to a FASTA file for use in ```Bifrost query```.

```python gfa_to_fasta.py graph.gfa```

The resulting FASTA file will take the same name as the graph, and can then be used in ```Bifrost query```. NOTE: exact k-mer matching should be used (```-e 1.0```).

```Bifrost query -e 1.0 -g graph.gfa -f graph.bfg_colors -q graph.fasta -o colours```

This will produce a TSV file named ```colours.tsv``` which can be used in ggCaller.

## Usage

To run ggCaller from the command line, pass arguments specifying input/output and parameters for traversal.
Arguments taken by ggcaller are:
- ```--graph``` Input GFA
- ```--colours``` Input TSV
- ```--kmer``` k-mer size used in ```Bifrost build```
- ```--path``` maximum path length in base-pairs (default 10000 bp)
- ```--orf``` minimum ORF length to return in base-pairs (default 90 bp)
- ```--out``` output file in fasta format (default 'calls.fasta')

In the below example, ggCaller is run against ```graph.gfa``` and ```colours.tsv```, querying a graph constructed using k=31. Maximum path length is set to 5000 bp, and minimum ORF length is set to 150 bp. A FASTA file called ```calls.fasta``` is then generated in the same directory.

```python ggcaller.py --graph graph.gfa --colours colours.tsv --kmer 31 --path 5000 --orf 150 --out calls.fasta```

## Interpreting output FASTA

The output ```calls.fasta``` contains ORF sequences, with headers which contain a unique identifier (```Gene_ID```), the strand of the ORF (```Strand```) and a colours array (```Colours```), which describes the presence/absence of the ORF in the source genomes in the same order as the matrix in ```colours.tsv```.

```
>[Gene_ID: 1] [Strand: +] [Colours: ['0', '0', '1', '1']]
TTGTCATTTTATTGGTTTATTTTCTTAGCTTTGTTAGAGAGACAGAACTTGAACGTTCTTC
>[Gene_ID: 2] [Strand: +] [Colours: ['1', '1', '1', '1']]
ATGGGATTAATGCTGTCTTTATGGGTGTTGGTGGCAGTTTTGATGTATTATCAGGACACATTAAACGAGCTCCATTATGGATGCAAAAATTGA
>[Gene_ID: 3] [Strand: +] [Colours: ['0', '0', '0', '1']]
TTGTTCAAAGGTGGTGTTACGATTTCAAGAACTCCTCTCAGTTCTGAGGACACGGTAATGATTGATGCGA
```

## Citation

If you use this code, please cite:

Bifrost: 
Holley G., Melsted, P. Bifrost – Highly parallel construction and indexing of colored and compacted de Bruijn graphs. bioRxiv 695338 (2019). doi: https://doi.org/10.1101/695338

pyGFA: 
https://github.com/AlgoLab/pygfa

Networkx: 
Hagberg AA, Schult DA, Swart PJ. Exploring Network Structure, Dynamics, and Function using NetworkX. In: Varoquaux G, Vaught T, Millman J, editors. Proceedings of the 7th Python in Science Conference. Pasadena, CA USA; 2008. p. 11–5




