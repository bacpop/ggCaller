Tutorial
==================================

Here we'll walk through a typical run of ggCaller, including both :ref:`Gene-calling` and :ref:`Querying`.

Installation and setup
----------------------

Follow the guide in :doc:`<installation>` for downloading and installing ggCaller.

Working Dataset
---------------

We'll use a dataset from `Bentley et al. (2006) <https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0020031>`_.
This dataset contains 91 sequences pneumococcal capsular polysaccharide synthetic (CPS) loci. These sequences are stucturally diverse,
but are only ~20,000 bp in length, so can be analysed quickly on a standard laptop or desktop.

Download the files from `here <>`_ and unzip::

    tar xvf Bentley_et_al_2006_CPS_sequences.tar.bz2

We will also provide out own custom annotation database for DIAMOND. These will be the manually curated protein sequences
from Bentley et al. Download from `here <>`_ and unzip::

    tar xvf Bentley_et_al_2006_CPS_protein_sequences.tar.bz2

Gene-calling
------------

First generate an input file for ggCaller. This must be a file containing paths (absolute recommended) to all sequences to be analysed.
We recommend running the below command within the unzipped to generate this file::

    cd Bentley_et_al_2006_CPS_sequences
    ls -d -1 $PWD/*.fa > input.txt
    cd ..

`input.txt` will now contain all `.fa` files in the directory `Bentley_et_al_2006_CPS_sequences`

Now we will run ggCaller specifying the below settings:

- Sensitive DIAMOND annotation using a custom database, and HMMER3 using the default database
- Pangenome-wide alignment using default MAFFT
- Saved intermediate datastructures, enabling sequence querying

To do this using 4 threads, run::

    ggcaller --refs fastas/Bentley_et_al_2006_CPS_sequences/input.txt --annotation ultrasensitive --diamonddb Bentley_et_al_2006_CPS_protein_sequences.faa --aligner def --alignment pan --save --out ggc_Bentley_et_al_CPS --threads 4

You will find the following files in the output directory ``ggc_Bentley_et_al_CPS``:

-