#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :

from setuptools import setup, find_packages
from codecs import open
from os import path
import os
import re
import io


def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


here = path.abspath(path.dirname(__file__))


with open(path.join(here, '..', 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='ggCaller',
    version=find_version("../ggCaller/__init__.py"),
    description='ggCaller: a bacterial gene caller for pangenome graphs',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/samhorsfield96/ggCaller',
    author='Sam Horsfield',
    author_email='shorsfield@ebi.ac.uk',
    license='MIT License',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
    ],
    python_requires='>=3.8.0',
    keywords='bacteria genomics pangenome genes annotation',
    packages=['ggCaller'],
    entry_points={'console_scripts': ['ggcaller = ggCaller.__main__:main']},
    test_suite="test"
)