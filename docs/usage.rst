Usage
==================================

ggCaller has two main modes: :ref:`Gene-calling` and :ref:`Querying`.

.. _gene-calling:

Gene-calling
-------------

Gene-calling predicts and annotates genes within a pangenome de Bruijn Graph (DBG), before
conducting orthologue clustering and pangenome analysis using Panaroo.

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

This assumes all sequences used to build the graph are assemblies.
If only some sequences are assemblies and the rest are reads, specify which files are references using ``--refs``::

    ggcaller --graph input.gfa --colours colours.bfg_colors --refs input1.txt

If all sequences are reads, specify ``--not-ref``::

    ggcaller --graph input.gfa --colours colours.bfg_colors --not-ref

Results from all commands above will be saved to a directory called ``ggCaller_output`` by default.
To change this, specify ``--out``. Note that ggCaller will overwrite results if an already existing directory is specified.

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

Quality control
^^^^^^^^^^^^^^^

ggCaller implements Panaroo to identify spurious clusters that are generated by assembly fragmentation and contamination.

Panaroo identifies spurious clusters as those with <2 edges in the gene graph. Spurious clusters are then removed based
on their population frequency, determined by three settings:

- ``strict``; remove spurious clusters with <5% frequency. Good for datasets >100 genomes where rare plasmids are not expected.
- ``moderate``; remove spurious clusters with <1% frequency (default). Good for datasets <=100 genomes where rare plasmids are not expected.
- ``sensitive``; do not remove clusters. Good for datasets where rare plasmids are expected.

For example, to run ggCaller in strict mode::

    ggcaller --refs input.txt --clean-mode strict

More information can be found `here <https://gtonkinhill.github.io/panaroo/#/gettingstarted/params>`_.

.. _querying:

Querying
--------

Querying maps a set of query DNA sequences to an annotated DBG, identifying genes that
the query overlaps with.

Saving datastructures
^^^^^^^^^^^^^^^^^^^^^

Annotate a DBG as before, adding the ``--save`` flag. This will write the intermediate datastructures
containing DBG coordinates of the predicted genes to a directory called ``ggc_data``.

.. important::
    We suggest using an annotation database, either the default
    ones provided or a custom one, as this will enable better
    functional analysis of your queries.

For example, run with sensitive annotation and save intermediate files::

    ggcaller --refs input.txt --annotation sensitive --save

Querying the DBG
^^^^^^^^^^^^^^^^^^^^^

Queries sequences can either be in multi-FASTA format, or in a single file with each sequence on its own line.

Provide paths to the DBG ``.gfa`` and ``.bfg_colors`` files, the ``ggc_data`` directory and query file::

    ggcaller --query queries.fasta --graph inputs.gfa --colours inputs.bfg_colors --data ggCaller_output/ggc_data

By default, mapped queries >=80% matching k-mers to a given colour will be returned. This can be changed using
``--query-id`` flag.

To return queries with 100% match::

    ggcaller --query queries.fasta --graph inputs.gfa --colours inputs.bfg_colors --data ggCaller_output/ggc_data --query-id 1.0

Interpreting results
^^^^^^^^^^^^^^^^^^^^^

Results will be output in ``matched_queries.fasta`` in the specified output directory. This is a multi-FASTA file describing
all annotated genes that overlap with the query sequences.

An example format is below::

    >Isolate10_9298 ggcID=10_9298 QUERY=TACTGCTAAACAAAGTCGCAAAGAAATCGAA;GAGTCTAGCTAGCTAATCGATAGCTACTAGA annotation=FUNCTION A;FUNCTION B;
    ATGTTAAATAAAGTCAAAACTAAAGCCTTAATTAGTGTCGGAGCAGTGGCTGCAACTAGCTAG

The header contains:

- Sample name and gene number (``Isolate10_9298``)
- ggCaller identifier (``ggcID`` field)
- Mapped query sequence(s) (``QUERY`` field) separated by semi-colons
- Annotation(s) (``annotation`` field) separated by semi-colons

Parallelisation
---------------

ggCaller is fully parallelised using OpenMP and python multiprocessing. By default ggCaller runs single-threaded.

To specify the number of threads::

    ggcaller --refs input.txt --threads 8
