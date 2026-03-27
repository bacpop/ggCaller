Usage
==================================

ggCaller has two main modes: :ref:`Gene-calling` and :ref:`Querying`.

.. _gene-calling:

Gene-calling
-------------

.. important::
    As of ggCaller v1.4, we highly recommend running ggCaller with
    ``--gene-finding-only`` enabled. This will not conduct annotation or
    COG clustering, and will instead generate a directory of GFF files
    which can be used in any modern clustering method, such as Panaroo.

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
    intermediate FMINDEX files in the same locations.

If not using Docker, generate the input file for ggCaller, navigate inside the directory containing the genomes, and run::

    ls -d -1 $PWD/*.fasta > input.txt

If using Docker, you must navigate to the directory containing the fasta files and run::

    ls -d -1 *.fasta > input.txt

This will generate a list of all the ``.fasta`` files in the directory. Change this extension as required.

.. important::
    All of the below commands can be run with docker installations, however they
    must be run as: ``docker run --rm -it -v $(pwd):/workdir samhorsfield96/ggcaller:latest ggcaller <commands>``.
    This command must be run within the same directory as the `.fasta` files and `input.txt`.
    All paths provided must be relative, as absolute paths will not work within the docker container.

DBG building with reads or assemblies is different, with k-mers that appear only once being removed from the graph.
Therefore it is important to specify whether ``input.txt`` contains reads or assemblies.

.. important::
    Assemblies with many Ns generate disjointed DBGs leading
    to underclustering. To ensure optimal performance, avoid
    using assemblies containing Ns.

To run ggCaller with just assemblies::

    ggcaller --refs input.txt

To run ggCaller with just reads::

    ggcaller --reads input.txt

To run ggCaller with reads and assemblies::

    ggcaller --refs input1.txt --reads input2.txt

.. important::
    We haven't extensively tested calling genes within
    read datasets yet. Exercise caution when interpreting
    results.

ggCaller can also be run on a pre-built Bifrost DBG and its associated colours file::

    ggcaller --graph input.gfa --colours colours.color.bfg

This assumes all sequences used to build the graph are assemblies.
If only some sequences are assemblies and the rest are reads, specify which files are references using ``--refs``::

    ggcaller --graph input.gfa --colours colours.color.bfg --refs input1.txt

If all sequences are reads, specify ``--not-ref``::

    ggcaller --graph input.gfa --colours colours.color.bfg --not-ref

You can also reduce runtime by skipping Panaroo with ``--gene-finding-only`` to generate just GFF and FASTA files, which can then be used in pangenome clustering methods, in the ``GFF`` subdirectory::

    ggcaller --refs input.txt --gene-finding-only

Results from all commands above will be saved to a directory called ``ggCaller_output`` by default.
To change this, specify ``--out <path>``. Note that ggCaller will overwrite results if an already existing directory is specified.

By default, ggCaller will generate:

- Predicted genes (nucleotide and amino-acid) in FASTA format
- GFF files for each input genome in a separate directory ``GFF``

Additionally, ggCaller generates some intermediate files:

- Two Bifrost files, a GFA file and COLORS file, with the same file path as ``input.txt``
- FMINDEX files for each of the sample FASTAs, placed in the ``Path_dir`` subdirectory.
- GRAPH INDEX files used for memory efficient gene calling and gene call updating, placed in the ``ORF_dir`` subdirectory.

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

For example, save intermediate files::

    ggcaller --refs input.txt  --save

Querying the DBG
^^^^^^^^^^^^^^^^^^^^^

Queries sequences can either be in multi-FASTA format, or in a single file with each sequence on its own line.

Provide paths to the DBG ``.gfa`` and ``.color.bfg`` files, the previous run directory::

    ggcaller --query queries.fasta --graph inputs.gfa --colours inputs.color.bfg --prev-run ggCaller_output

By default, mapped queries >=80% matching k-mers to a given colour will be returned. This can be changed using
``--query-id`` flag.

To return queries with 100% match::

    ggcaller --query queries.fasta --graph inputs.gfa --colours inputs.color.bfg --prev-run ggCaller_output --query-id 1.0

.. _Interpreting results:

Interpreting results
^^^^^^^^^^^^^^^^^^^^^

Results will be output in ``matched_queries.fasta`` in the specified output directory. This is a multi-FASTA file describing
all annotated genes that overlap with the query sequences.

An example format is below::

    >Isolate10_9298 ggcID=10_9298 QUERY=Query_A;Query_B;
    ATGTTAAATAAAGTCAAAACTAAAGCCTTAATTAGTGTCGGAGCAGTGGCTGCAACTAGCTAG

The header contains:

- Sample name and gene number (``Isolate10_9298``)
- ggCaller identifier (``ggcID`` field)
- Mapped query sequences or IDs (``QUERY`` field) separated by semi-colons. These will be fasta IDs if ``queries`` file is a FASTA, otherwise DNA sequence.

Iterative gene calling
----------------------

After an initial run of ggCaller, you can call genes in new genomes, using the original information from the initial gene calls. 

This is designed to be used after a run with ``--gene-finding-only``, as it does not use information from Panaroo::

    ggcaller --refs input1.txt --gene-finding-only --out run1
    ggcaller --refs input2.txt --gene-finding-only --out run2 --prev-run run1

Results can be placed in a new directory, or directed to the original directory. If repeated updates are likely, use a single directory::

    ggcaller --refs input1.txt --gene-finding-only --out all_runs
    ggcaller --refs input2.txt --gene-finding-only --out all_runs --prev-run all_runs
    ggcaller --refs input3.txt --gene-finding-only --out all_runs --prev-run all_runs

Parallelisation
---------------

ggCaller is fully parallelised using OpenMP and python multiprocessing. By default ggCaller runs single-threaded.

To specify the number of threads::

    ggcaller --refs input.txt --threads 8
