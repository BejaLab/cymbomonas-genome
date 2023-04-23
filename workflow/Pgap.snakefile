
from Bio import SeqIO
import yaml
from os.path import basename

fasta = SeqIO.parse("contamination/bacteria.fna", "fasta")
records = SeqIO.to_dict(fasta)

rule all:
    input:
        expand("analysis/contamination/pgap/{id}/annot.tbl", id = records.keys())

rule fetch_seq:
    input:
        "contamination/bacteria.fna"
    output:
        "analysis/contamination/pgap/{id}.fna"
    conda:
        "envs/tools.yaml"
    shell:
        "seqkit grep -p {wildcards.id} -o {output} {input}"

rule copy_submol:
    input:
        "workflow/metadata/submol.yaml"
    output:
        "analysis/contamination/pgap/{id}_submol.yaml"
    shell:
        "cp {input} {output}"

rule make_yaml:
    input:
        fasta  = "analysis/contamination/pgap/{id}.fna",
        submol = "analysis/contamination/pgap/{id}_submol.yaml"
    output:
        "analysis/contamination/pgap/{id}.yaml"
    wildcard_constraints:
        id = '.*(?<!_submol)'
    run:
        sets = dict(
            fasta  = { "class": "File", "location": basename(input.fasta)  },
            submol = { "class": "File", "location": basename(input.submol) }
        )
        with open(output[0], 'w') as fd:
            yaml.dump(sets, fd)

rule pgap:
    input:
        yaml = "analysis/contamination/pgap/{id}.yaml",
        submol = "analysis/contamination/pgap/{id}_submol.yaml"
    output:
        gbk = "analysis/contamination/pgap/{id}/annot.gbk"
    params:
        dirname = directory("analysis/contamination/pgap/{id}/")
    threads:
        8
    shell:
        """
        rm -r {params.dirname}
        pgap.py --debug --ignore-all-errors --report-usage-true --cpu {threads} -o {params.dirname} {input.yaml}
        """

rule gbf2tbl:
    input:
        script = "workflow/tools/gbf2tbl.pl",
        pgap   = "analysis/contamination/pgap/{id}/annot.gbk"
    output:
        "analysis/contamination/pgap/{id}/annot.tbl"
    conda:
        "envs/perl.yaml"
    shell:
        "perl {input.script} {input.pgap}"
