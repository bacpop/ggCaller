Quickstart
==================================

Installation
------------

The easiest way to install is through conda, which will also install the
dependencies::

    conda install ggcaller

.. important::
    ggCaller requires python3.8 to run
    (which on many default Linux installations is
    run using ``python3`` rather than ``python``).

Preparing the data
------------------

Place all of you samples to analyse in the same directory. Then navigate inside and run::

    ls -d -1 $PWD/*.fasta > input.txt

Running ggCaller
------------------

To run ggCaller with just assemblies::

    ggcaller --refs input.txt

To run ggCaller with just reads::

    ggcaller --reads input.txt

Results will be saved to the directory ``ggCaller_output`` by default. To change this, specify ``--out <path>``.
