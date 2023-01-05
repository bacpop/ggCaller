Usage
==================================

ggCaller has two main modes; gene-calling and querying.

Gene-calling predicts and annotates genes within a pangenome de Bruijn Graph (DBG), before
conducting orthologue clustering and pangenome analysis using Panaroo.

Querying maps a set of query DNA sequences to an annotated DBG, identifying genes that
the query overlaps with.

Gene calling
------------

Predicting genes
^^^^^^^^^^^^^^^^

To generate an input for ggCaller, create a directory containing of all the sequences you wish to analyses.
We recommend placing all samples of the same type in a single directory; place read and assembly files in
separate directories.

.. important::
    Ensure you have write access to the directories where
    the FASTA/FASTQ files are saved, as ggCaller saves
    intermediate files in the same locations.

To generate the input file for ggCaller, navigate inside the directory containing the genomes, and run::

    ls -d -1 $PWD/*.fasta > input.txt

This will generate a list of all the ``.fasta`` files in the directory. Change this extension as required.

DBG building with reads or assemblies is different, with k-mers that appear only once being removed from the graph.
Therefore it is important to specify whether ``input.txt`` contains reads or assemblies.

To run ggCaller with just assemblies::

    ggcaller --refs input.txt

To run ggCaller with just reads::

    ggcaller --reads input.txt

To run ggCaller with reads and assemblies::

    ggcaller --refs input1.txt --reads input2.txt

ggCaller can also be run on a pre-built Bifrost DBG and its associated colours file::

    ggcaller --graph input.gfa --colours colours.bfg_colors

This will save results to a directory called ``ggCaller_output`` by default. To change this, specify ``--out``.
Note that ggCaller will overwrite results if an already existing directory is specified.

By default, ggCaller will generate:
- Predicted genes (nucleotide and amino-acid) in FASTA format
- Gene presence/absence matrix in CSV and RTAB formats
- Pre/post Panaroo quality control gene graphs in GML format
- Structural variant presence/absence in RTAB format
- Summary graph: gene frequency, cluster size and rarefaction curve
- Roary-style gene frequency statistics
- A pangenome reference FASTA, containing all cluster centroids
- A gene presence/absence neighbour joining tree in NWK format

Annotating genes
^^^^^^^^^^^^^^^^

ggCaller comes with two default databases for functional annotation of genes.
- Bacterial and Viral databases from `Uniprot <https://www.uniprot.org/>`_, used by `DIAMOND <https://github.com/bbuchfink/diamond>`_
- HMM profiles from `Prokka <https://github.com/tseemann/prokka>`_, used by `HMMER3 <https://github.com/EddyRivasLab/hmmer>`_

.. important::
    Ensure you are connected to the internet
    when first running ggCaller as these databases
    are downloaded automatically. Subsequent runs
    can be conducted offline.

There are three sensitivity levels for annotation:
- ``fast``: only DIAMOND  in fast mode
- ``sensitive``: only DIAMOND in sensitive mode
- ``ultrasensitive``: HMMER3 and DIAMOND in sensitive mode

For example, to run DIAMOND only in fast mode, run::

    ggcaller --refs input.txt --annotation fast

By default these commands will annotate using DIAMOND with the ``Bacteria`` uniprot database.
To change this to the ``Viruses`` database, run::

    ggcaller --refs input.txt --annotation fast --diamonddb Viruses

Custom databases can also be specified for both DIAMOND using ``--diamonddb`` and HMMER3 using ``--hmmdb``.
DIAMOND databases must be amino-acid FASTA files. HMMER3 databases must be HMM-profile ``.HAMAP`` files built using
``hmmbuild`` which is part of the HMMER3 package.

To run with custom DIAMOND and HMMER3 databases::

    ggcaller --refs input.txt --annotation ultrasensitive --diamonddb annotation.fasta --hmmdb annotation.HAMAP

Annotation is not on by default. If annotation is specified, ggCaller will additionally generate:
- GFF files for each input genome in a separate directory ``GFF``
- Annotations will be added to gene call FASTA files

Aligning genes
^^^^^^^^^^^^^^

ggCaller also supports generation of within-cluster and core genome alignments using `MAFFT <https://github.com/GSLBiotech/mafft>`_.

There are two alignment algorithms implemented:
- ``def`` or default, which uses the standard MAFFT multiple sequence alignment algorithm. This is faster when aligning <=500 sequences in a cluster.
- ``ref`` or reference, which uses reference-guided alignment. This is faster when aligning >500 sequences in a cluster.

There are also two modes for alignment:
- ``core`` aligns genes only within core clusters, and generates a concatenated core genome alignment.
- ``pan`` aligns genes within all clusters (pangenome alignment), as well as generating a concatenated core genome alignment.

To generate a core genome alignment  using default MAFFT, run::

    ggcaller --refs input.txt --aligner def --alignment core

To generate a pangenome alignment using reference-guided MAFFT, run::

    ggcaller --refs input.txt --aligner ref --alignment pan

To change the frequency of genes deemed to be core, use `--core-threshold` (default = 0.95, or 95% frequency).
For example, only include genes found at 100% frequency::

    ggcaller --refs input.txt --aligner def --alignment core --core-threshold 1.0

Alignment is off by default. If specified, ggCaller will additionally generate:
- Core genome alignment in FASTA format
- Core genome Neighbour-joining tree in NWK format
- Per-cluster alignment files in FASTA format in a separate directory ``aligned_gene_sequences``
- Per-cluster VCF file generated by `SNP-SITES <https://github.com/sanger-pathogens/snp-sites>`_ in separate directory ``VCF``