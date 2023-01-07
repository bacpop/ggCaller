Installation
============

ggCaller is available on Linux and MacOS. If you are running Windows 10/11, Linux can be installed via the Windows Subsystem for Linux (`WSL <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_).

The easiest way to install is through conda, which will also install the
dependencies::

    conda install ggcaller

.. important::
    ggCaller requires python3.8 to run
    (which on many default Linux installations is
    run using ``python3`` rather than ``python``).

Installing with conda (recommended)
-----------------------------------
If you do not have ``conda`` you can install it through
`miniconda <https://conda.io/miniconda.html>`_ and then add the necessary
channels::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

Then run::

    conda install ggcaller

Clone the code
--------------
You can also clone the github to run the latest version, which is executed by::

    git clone --recursive https://github.com/samhorsfield96/ggCaller && cd ggCaller
    python setup.py install

You will need to install the dependencies yourself (you can still use
conda or pip for this purpose). See ``environment.yml``. In addition,
a C++17 compiler (e.g. gcc >=7.3) is required.

Test installation
-----------------

After any of the above steps, check correct setup by running ``ggcaller --help``.