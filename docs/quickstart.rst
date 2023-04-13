Quickstart
==================================

.. important::
    We are aware of issues installing from conda version at the moment.
    We recommend installing from source at this time.

Installation
------------

The easiest way to get up and running is using . First, download and switch to the ggCaller repository::

    git clone --recursive https://github.com/samhorsfield96/ggCaller && cd ggCaller

Then, build with Docker. This should take between 5-10 minutes to fully install.::

	docker build -t ggc_env:latest -f docker/Dockerfile .


Preparing the data
------------------

Place all of you samples to analyse in the same directory. Then navigate inside and run::

    ls -d -1 $PWD/*.fasta > input.txt

If using Docker, instead navigate to the directory containing the fasta files and run the below command, to ensure paths files are relative. The docker version will not work with absolute paths.::

    ls -d -1 *.fasta > input.txt

Running ggCaller
------------------

To run ggCaller with just assemblies::

    ggcaller --refs input.txt --out output_path

To run ggCaller with just reads::

    ggcaller --reads input.txt --out output_path

If using Docker, run with the below command. You must ensure all paths are relative, including in ``input.txt``::

	docker run --rm -it -v $(pwd):/workdir ggc_env:latest ggcaller --refs input.txt --out output_path

.. important::
    We haven't extensively tested calling genes within
    read datasets yet. Exercise caution when interpreting
    results.

Results will be saved to the directory ``ggCaller_output`` by default. To change this, specify ``--out <path>``.
