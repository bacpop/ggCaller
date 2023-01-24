Installation
============

ggCaller is available on Linux. If you are running Windows 10/11, Linux can be installed via the Windows Subsystem for Linux (`WSL <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_).

We plan to get a MacOS version up and running in the future.

.. important::
    ggCaller requires python3.8 to run
    (which on many default Linux installations is
    run using ``python3`` rather than ``python``).

Installing with conda
-----------------------------------

.. important::
    We are aware of issues installing from conda at the moment.
    We recommend installing from source at this time.

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

For example, using conda (creates ```ggc_env``` environment)::

    conda env create -f environment_linux.yml
    conda activate ggc_env

Then clone the code and install::

    git clone --recursive https://github.com/samhorsfield96/ggCaller && cd ggCaller
    python setup.py install

Test installation
-----------------

After any of the above steps, check correct setup by running ``ggcaller --help``.