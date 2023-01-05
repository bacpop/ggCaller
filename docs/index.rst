ggCaller documentation
==================================
ggCaller is a novel bacterial gene annotation and pangenome analysis tool, designed to enable fast, accurate analysis of large single-species genome datasets.

ggCaller traverses de Bruijn graphs (DBGs) built by `Bifrost <https://github.com/pmelsted/bifrost>`_,
using temporal convolutional networks from `Balrog <https://github.com/salzberg-lab/Balrog>`_ for gene filtering
and `Panaroo <https://github.com/gtonkinhill/panaroo>`_ for pangenome analysis and quality control.

.. image:: images/ggCaller_logo.png
   :alt:  ggCaller (Graph Gene Caller)
   :align: center

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   self
   installation.rst

Why ggCaller?
-------------
ggCaller uses population-frequency information at several stages of gene annotation and pangenome analysis. This has several benefits:

- Consistent identification of start and stop codons across orthologs, improving clustering accuracy.
- Reduced gene-annotation sensitivity to assembly fragmentation.
- Reduced runtime verses existing gene-annotation and pangenome analysis workflows.
- One-line command from fasta -> gene annotations, gene frequency matrices, clusters of orthologous genes (COGs), core genome/pangenome alignments, phylogenetic trees, small/structural variants and more!
- Annotated DBG-querying for functional PanGenome-Wide Association Studies (PGWAS), compatible with results from `Pyseer <https://github.com/mgalardini/pyseer>`_.

Contents
---------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
