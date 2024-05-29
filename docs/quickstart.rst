Quickstart
==================================

.. important::
    We are aware of issues installing from conda version at the moment.
    We recommend installing via Docker at this time.

Installation
------------

The easiest way to get up and running is using Docker. To get up and running, pull the latest image::

    docker pull samhorsfield96/ggcaller:latest

Preparing the data
------------------

Place all of your samples to be analysed in the same directory. Then navigate inside and run::

    ls -d -1 $PWD/*.fasta > input.txt

If using Docker, instead navigate to the directory containing the fasta files and run the below command, to ensure file paths are relative (the docker version will not work with absolute paths)::

    ls -d -1 *.fasta > input.txt

Running ggCaller
------------------

.. important::
    Assemblies with many Ns generate disjointed DBGs leading
    to underclustering. To ensure optimal performance, avoid
    using assemblies containing Ns.

To run ggCaller with just assemblies::

    ggcaller --refs input.txt --out output_path

To run ggCaller with just reads::

    ggcaller --reads input.txt --out output_path

If using Docker, run with the below command. You must ensure all paths are relative, including in ``input.txt``::

	docker run --rm -it -v $(pwd):/workdir samhorsfield96/ggcaller:latest ggcaller --refs input.txt --out output_path

.. important::
    We haven't extensively tested calling genes within
    read datasets yet. Exercise caution when interpreting
    results.

Results will be saved to the directory ``ggCaller_output`` by default. To change this, specify ``--out <path>``.
