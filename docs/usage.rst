Usage
==================================

ggCaller has two main modes; gene-calling and querying.

Gene-calling predicts and annotates genes within a pangenome de Bruijn Graph (DBG), before
conducting orthologue clustering and pangenome analysis using Panaroo.

Querying maps a set of query DNA sequences to an annotated DBG, identifying genes that
the query overlaps with.

Gene calling
------------

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

This will save results to a directory called ``ggCaller_output`` by default. To change this, specify ``--out``.
Note that ggCaller will overwrite results if an already existing directory is specified.
