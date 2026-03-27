Advanced
==================================

For advanced users, ggCaller has a number of parameters for altering gene prediction, annotation and quality control.

Input/output
^^^^^^^^^^^^

- ``--kmer``: value of k used to build Bifrost DBG (Default and max value = 31).
- ``--all-seq-in-graph``: Output gene graph GML file with all DNA and amino acid sequences. Off by default due to large file size.
- ``--balrog-db``: Path to an existing download of the balrog annotation database. If this does not exist, downloaded and placed in path for future use.
- ``--gene-finding-only``: Only run ggCaller gene-finding and generate a gff compatible with other clustering tools.

Traversal and gene-calling cut-off settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``--max-path-length``: Maximum path length traversed during ORF finding (bp) (Default = 20000)
- ``--min-orf-length``: Minimum ORF length to return (bp) (Default = 90)
- ``--score-tolerance``: Probability threshold for shorter alternative start sites based on average stop codon frequency (Default = 0.2)
- ``--max-ORF-overlap``: Maximum overlap allowed between two ORFs (bp) (Default = 60)
- ``--min-path-score``: Minimum total BALROG score for a maximum tiling path of ORFs to be returned (Default = 100)
- ``--min-orf-score``: Minimum individual Balrog score for an ORF to be returned (Default = 100)
- ``--max-orf-orf-distance``: Maximum distance between two ORFs to be connected (bp) (Default = 10000)

Avoid/include algorithms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``--no-filter``: Do not filter ORF calls using Balrog, will return all ORF calls (Default = False)
- ``--no-write-idx``: Do not write FMIndexes to file (Default = False)
- ``--no-write-graph``: Do not write Bifrost GFA and colours to file (Default = False)
- ``--repeat``: Enable traversal of nodes multiple times, only applicable when DBG built from reads (Default = False)

Misc. options
^^^^^^^^^^^^^^^^^^

- ``--quiet``: Suppress additional output to console (Default = False)
- ``--version``: Show program's version number and exit (Default = False)