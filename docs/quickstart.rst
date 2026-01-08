Quickstart
==================================

Installation
------------

The easiest way to get up and running is using Docker. To get up and running, pull the latest image::

    docker pull samhorsfield96/ggcaller:latest

Preparing the data
------------------

Place all of your samples to be analysed in the same directory. Then navigate inside and run::

    ls -d -1 $PWD/*.fasta > input.txt

If using Docker, instead navigate to the directory containing the fasta files and run the below command, to ensure file paths are relative (the docker version will not work with absolute paths)::

    ls -d -1 *.fasta > input_docker.txt

Then, append the prefix ``/data/`` to each line to enable ggCaller to find the files::

    sed -i -e 's|^|/data/|' input_docker.txt

Running ggCaller
------------------

.. important::
    Assemblies with many Ns generate disjointed DBGs leading
    to underclustering. To ensure optimal performance, avoid
    using assemblies containing Ns.

To run ggCaller with just assemblies::

    ggcaller --refs input.txt --out output_path --gene-finding-only

.. important::
    As of ggCaller v1.4, we highly recommend running ggCaller with
    ``--gene-finding-only`` enabled. This will not conduct annotation or
    COG clustering, and will instead generate a directory of GFF files
    which can be used in any modern clustering method, such as Panaroo.

To run ggCaller with just reads::

    ggcaller --reads input.txt --out output_path --gene-finding-only

If using Docker, run with the below command, ensuring you keep ``--balrog-db /app/ggc_db`` and ``/workdir`` paths as specified below. Replace ``path to files`` with the absolute path to the directory of files in ``input_docker.txt``::

	docker run --rm -it -v $(pwd):/workdir -v <path to files>:/data samhorsfield96/ggcaller:latest ggcaller --balrog-db /app/ggc_db --refs /workdir/input_docker.txt --out /workdir/output_path --gene-finding-only

.. important::
    We haven't extensively tested calling genes within
    read datasets yet. Exercise caution when interpreting
    results.

Results will be saved to the directory ``ggCaller_output`` by default. To change this, specify ``--out <path>``.
