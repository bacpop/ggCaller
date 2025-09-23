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
- ``--no-clustering``: Do not cluster ORFs (Default = False)
- ``--no-refind``: Do not refind missed genes (Default = False)

Gene clustering options
^^^^^^^^^^^^^^^^^^^^^^^^

- ``--identity-cutoff``: Minimum identity at amino acid level between two ORFs for lowest-level clustering (Default = 0.98)
- ``--len-diff-cutoff``: Minimum ratio of length between two ORFs for lowest-level clustering (Default = 0.98)
- ``--family-threshold``: Gene family sequence identity threshold (default=0.7)
- ``--merge-paralogs``: Don't split paralogs during Panaroo quality control (Default = False)

Annotation options
^^^^^^^^^^^^^^^^^^^

- ``--evalue``: Maximum e-value to return for DIAMOND and HMMER searches during annotation (Default = 0.001)
- ``--truncation-threshold``: Sequences in a cluster less than `centroid length * truncation-threshold` will be annotated as 'potential pseudogene' (Default = 0.8)

Gene-refinding options
^^^^^^^^^^^^^^^^^^^^^^^

- ``--search-radius``: The distance (bp) surrounding the neighbour of an accessory gene in which to search for it (Default = 5000)
- ``--refind-prop-match``: The proportion of an accessory gene's length that must be found in order to consider it a match (Default = 0.2)

Gene graph correction stringency options (determined by clean-mode)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``--min-trailing-support``: Minimum cluster size to keep a gene called at the end of a contig.
- ``--trailing-recursive``: Number of times to perform recursive trimming of low support nodes near the end of contigs
- ``--edge-support-threshold``: Minimum support required to keep an edge that has been flagged as a possible mis-assembly
- ``--length-outlier-support-proportion``: Proportion of genomes supporting a spurious long gene (>1.5x outside the IQR of cluster)
- ``--min-edge-support-sv``: Minimum edge support required to call structural variants in the presence/absence sv file
- ``--no-clean-edges``: Turn off edge filtering in the final output graph

Alignment options
^^^^^^^^^^^^^^^^^^

- ``--no-variants``: Do not call variants using SNP-sites after alignment (Default = False)
- ``--ignore-pseduogenes``: Ignore ORFs annotated as 'potential pseudogenes' in alignments (Default = False)

Misc. options
^^^^^^^^^^^^^^^^^^

- ``--quiet``: Suppress additional output to console (Default = False)
- ``--version``: Show program's version number and exit (Default = False)