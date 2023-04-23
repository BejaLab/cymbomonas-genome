
import Bio

assembly_file = str(snakemake.input['assembly'])
red_bed_file  = str(snakemake.input['red_bed'])
homologs_file = str(snakemake.input['homologs'])
viruses_file  = str(snakemake.input['viruses'])
contamin_file = str(snakemake.input['contamin'])

