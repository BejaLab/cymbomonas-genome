
import yaml
from os.path import basename

ids ,= glob_wildcards("bacteria/{id}.fna")

rule all:
    input:
        expand("analysis/contamination/pgap/{id}/annot.tbl", id = ids)

rule copy_seq:
    input:
        "bacteria/{id}.fna"
    output:
        "analysis/contamination/pgap/{id}.fna"
    shell:
        "cp {input} {output}"

rule copy_submol:
    input:
        "bacteria/{id}_submol.yaml"
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
        rm -fr {params.dirname}
        pgap.py --no-internet --no-self-update --debug --ignore-all-errors --report-usage-true --cpu {threads} -o {params.dirname} {input.yaml}
        """

rule gbf2tbl:
    input:
        script = "workflow/lib/gbf2tbl.pl",
        pgap   = "analysis/contamination/pgap/{id}/annot.gbk"
    output:
        "analysis/contamination/pgap/{id}/annot.tbl"
    conda:
        "envs/perl.yaml"
    shell:
        "perl {input.script} {input.pgap}"
