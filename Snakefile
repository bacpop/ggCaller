import os

WORKDIR = os.getenv("PANWORKDIR", "./")

rule diamond_makedb:
    input:
       {ref_fasta}
    output:
        WORKDIR/{reference}.db
    shell:
        "diamond makedb --in {input} -d {output}"

rule diamond_search_90perc:
    input:
        reference = {reference}.db,
        query = WORKDIR/annotation/{samples}.fasta
    output:
        WORKDIR/annotation/{samples}.tsv
    shell:
        "diamond blastp --fast -d {input.reference} -q {input.query} -o {output}"

rule diamond_search_70perc:
    input:
        reference = {reference}.db,
        query = WORKDIR/annotation/{samples}.fasta
    output:
        WORKDIR/annotation/{samples}.tsv
    shell:
        "diamond blastp --mid-sensitive -d {input.reference} -q {input.query} -o {output}"

rule diamond_search_0perc:
    input:
        reference = {reference}.db,
        query = WORKDIR/annotation/{samples}.fasta
    output:
        WORKDIR/annotation/{samples}.tsv
    shell:
        "diamond blastp --ultra-sensitive -d {input.reference} -q {input.query} -o {output}"