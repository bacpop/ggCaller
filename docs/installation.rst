Installation
============

ggCaller is available on Linux. If you are running Windows 10/11, Linux can be installed via the Windows Subsystem for Linux (`WSL <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_). If running via Docker, ensure you install WSL2.

We plan to get a MacOS version up and running in the future.

.. important::
    ggCaller requires python3.9 to run
    (which on many default Linux installations is
    run using ``python3`` rather than ``python``).

Installing with Docker
-----------------------------------

First, install `Docker <https://docs.docker.com/get-docker/>`_ for your OS. If running with WSL2, you should still download Docker Desktop for Windows.

To use the latest image, run::

    docker pull samhorsfield96/ggcaller:latest

To run ggCaller from the Docker Hub image, run::

	cd test && docker run --rm -it -v $(pwd):/workdir -v $(pwd):/data samhorsfield96/ggcaller:latest ggcaller --balrog-db /app/ggc_db --refs /workdir/pneumo_CL_group2_docker.txt --out /workdir/ggc_out

You can also build the image yourself. First download and switch to the ggCaller repository::

    git clone --recursive https://github.com/samhorsfield96/ggCaller && cd ggCaller

Finally, build with Docker. This should take between 5-10 minutes to fully install.::

	docker build -t ggc_env:latest -f docker/Dockerfile .

To run ggCaller from a local Docker build, run::

	cd test && docker run --rm -it -v $(pwd):/workdir -v $(pwd):/data ggc_env:latest ggcaller --balrog-db /app/ggc_db --refs /workdir/pneumo_CL_group2_docker.txt --out /workdir/ggc_out

Please ensure you keep ``--balrog-db /app/ggc_db`` and ``/workdir`` paths as specified above.

Installing with singularity
-----------------------------------

If you encounter permissions issues using Docker, you can download the singularity image from `Zenodo <https://zenodo.org/record/7870950>`_

Once downloaded, set up the singularity container using::

    singularity shell --writable <singulatiry image>.sif

Once loaded, add the conda bin directory to your path variable and run ggCaller as normal::

    PATH=$PATH:/opt/conda/bin
    ggcaller --refs input.txt --out output_path

Installing with conda
-----------------------------------

Installing with conda is the easiest way to get ggCaller up and running, and will install all dependencies.

If you do not have ``conda`` you can install it through
`miniconda <https://conda.io/miniconda.html>`_ and then add the necessary
channels::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

Then run::

    conda install ggcaller

Installing from source
----------------------

You can also clone the github code to run the latest version.

You will need to install the dependencies yourself (you can still use
conda or pip for this purpose). See ``environment_linux.yml`` or ``environment_macOS.yml``.
In addition, a C++17 compiler (e.g. gcc >=7.3) is required.

We highly recommend using mamba over conda due to the large number of dependencies, making mamba significantly faster.

To install dependencies (creates ```ggc_env``` environment)::

    mamba env create -f environment_linux.yml
    mamba activate ggc_env

Then clone the code and install::

    git clone --recursive https://github.com/samhorsfield96/ggCaller && cd ggCaller
    python setup.py install

Test installation
-----------------

After any of the above steps, check correct setup by running ``ggcaller --help``.