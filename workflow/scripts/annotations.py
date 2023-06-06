
from Bio import SeqIO
from collections import defaultdict
import re

tbl_file = snakemake.input['tbl']
trna_file = snakemake.input['trna']
rrna_file = snakemake.input['rrna']
fsa_file = snakemake.input['fsa']
contam_file = snakemake.input['contamination']

contam_threshold = snakemake.params['contam_threshold']

trna_threshold = snakemake.params['trna_threshold']
trna_software = snakemake.params['trna_software']

rrna_threshold = snakemake.params['rrna_threshold']
rrna_software = snakemake.params['rrna_software']
rrna_database = snakemake.params['rrna_database']
species = snakemake.params['species']
strain = snakemake.params['strain']

tbl_out_file = snakemake.output['tbl']
fsa_out_file = snakemake.output['fsa']

def locus_tag_gen(species):
    GEN, SPEC, *rest = species.upper().split()
    lt_prefix = GEN[0:3] + SPEC[0:3] + '_'
    locus_num = 0
    while True:
        locus_num += 1
        yield lt_prefix + str(locus_num)

def read_sec(fh):
    re1 = re.compile('([^ ]+) \((\d+)-(\d+)\)')
    re2 = re.compile('Type: (\w+)\tAnticodon: (\w+) at \d+-\d+ \((\d+)-(\d+)\)\tScore: ([\d.]+)')
    re3 = re.compile('intron: \d+-\d+ \((\d+)-(\d+)\)')
    data = ''
    for line in fh:
        if not line.strip():
            locus_tag, start, end = re1.search(data).groups()
            aa, anticodon, ac_start, ac_end, score = re2.search(data).groups()
            is_pseudo = 'pseudogene' in data
            coords = [ (int(start), int(end)) ]
            forward = coords[0][0] < coords[0][1]
            if search := re3.search(data):
                i_start = int(search.group(1))
                i_end   = int(search.group(2))
                if forward:
                    coords = [ (coords[0][0], i_start - 1), (i_end + 1, coords[0][1]) ]
                else:
                    coords = [ (coords[0][0], i_start + 1), (i_end - 1, coords[0][1]) ]
            yield locus_tag, coords, aa, anticodon, int(ac_start), int(ac_end), float(score), is_pseudo
            data = ''
        else:
            data += line
    if data:
        yield data

def wrap_coords(feature, coords):
    start, end = coords[0]
    lines = [ "\t".join([ str(start), str(end), feature ]) ]
    for start, end in coords[1:]:
        lines.append("\t".join([ str(start), str(end) ]))
    return lines

def wrap_desc(desc):
    lines = []
    for key, vals in desc.items():
        if not isinstance(vals, list):
            vals = [ vals ]
        for val in vals:
            lines.append("\t".join([ "", "", "", key, str(val) ]))
    return lines

locus_tags = locus_tag_gen(species)

ncRNA = defaultdict(list)
with open(trna_file) as fh:
    for locus, trna_coords, aa, anticodon, ac_start, ac_end, score, is_pseudo in read_sec(fh):
        if score > trna_threshold and not is_pseudo:
            locus_tag = next(locus_tags)
            seqname = locus.rsplit(".", maxsplit = 1)[0]
            if ac_start < ac_end:
                ac_pos = '{start}..{end}'.format(start = ac_start, end = ac_end)
            else:
                ac_pos = 'complement({end}..{start})'.format(start = ac_start, end = ac_end)

            gene_desc = { "locus_tag": locus_tag }

            trna_desc = {
                "locus_tag": locus_tag,
                "inference": "COORDINATES: profile:{software}".format(software = trna_software)
            }
            if is_pseudo:
                gene_desc["pseudo"] = True
                trna_desc["pseudo"] = True
            else:
                if aa == "Undet":
                    trna_desc["product"] = "Undetermined tRNA"
                elif aa == "Sup":
                    trna_desc["product"] = "Suppressor tRNA"
                else:
                    trna_desc["product"] = "tRNA-" + aa
                    trna_desc["anticodon"] = "(pos:{pos},aa:{aa},seq:{anticodon})".format(pos = ac_pos, aa = aa, anticodon = anticodon.lower())

            feature_type = "tRNA" if aa not in [ "Undet", "Sup" ] else "ncRNA"

            lines = []
            lines += wrap_coords("gene", [ (trna_coords[0][0], trna_coords[-1][1]) ])
            lines += wrap_desc(gene_desc)
            if aa in [ "Undet", "Sup" ]:
                lines += wrap_coords("ncRNA", trna_coords)
                lines += wrap_desc({ "ncRNA_class": "other" })
            else:
                lines += wrap_coords("tRNA", trna_coords)
            lines += wrap_desc(trna_desc)
            ncRNA[seqname] += lines

products = {
    "SSU_rRNA_bacteria": ("16S ribosomal RNA", 1000),
    "SSU_rRNA_eukarya": ("18S ribosomal RNA", 1000),
    "LSU_rRNA_bacteria": ("23S ribosomal RNA", 2000),
    "LSU_rRNA_eukarya": ("28S ribosomal RNA", 3300),
    "5S_rRNA": ("5S ribosomal RNA", 90),
    "5_8S_rRNA": ("5.8S ribosomal RNA", 130)
}

with open(rrna_file) as fh:
    for line in fh:
        if not line.startswith('#'):
            idx, target_name, target_accession, query_name, query_accession, clan, mdl, mdl_from, mdl_to, seq_from, seq_to, strand, trunc, pass_, gc, bias, score, e_value, *rest = line.split()
            if target_name in products and float(e_value) < rrna_threshold:
                lines = []
                product, min_len = products[target_name]
                if abs(int(seq_to) - int(seq_from)) + 1 < min_len:
                    rrna_desc = {
                        "note": "Fragment of {product}".format(product = product),
                        "inference": [
                            "COORDINATES: profile:{software}".format(software = rrna_software),
                            "COORDINATES: nucleotide motif:{database}:{acc}".format(database = rrna_database, acc = target_accession)
                        ]
                    }
                    lines += wrap_coords("misc_feature", [ (seq_from, seq_to) ])
                    lines += wrap_desc(rrna_desc)
                else:
                    locus_tag = next(locus_tags)
                    gene_desc = {
                        "locus_tag": locus_tag
                    }
                    rrna_desc = {
                        "locus_tag": locus_tag,
                        "product": product,
                        "inference": [
                            "COORDINATES: profile:{software}".format(software = rrna_software),
                            "COORDINATES: nucleotide motif:{database}:{acc}".format(database = rrna_database, acc = target_accession)
                        ],
                        "db_xref": "RFAM:{acc}".format(acc = target_accession)
                    }
                    lines += wrap_coords("gene", [ (seq_from, seq_to) ])
                    lines += wrap_desc(gene_desc)
                    lines += wrap_coords("rRNA", [ (seq_from, seq_to) ])
                    lines += wrap_desc(rrna_desc)
                ncRNA[query_name] += lines

contamination = defaultdict(int)
with open(contam_file) as fh:
    for line in fh:
        seqname, offset, end = line.split()
        contamination[seqname] += int(end) - int(offset)

skip_sequences = {}
with open(fsa_file) as fsa:
    with open(fsa_out_file, 'w') as out:
        for record in SeqIO.parse(fsa, 'fasta'):
            skip_sequence = contamination[record.id] / len(record.seq) > contam_threshold
            if skip_sequence:
                skip_sequences[record.id] = True
            else:
                record.description = '[organism={species}] [strain={strain}]'.format(species = species, strain = strain)
                SeqIO.write(record, out, 'fasta')

def get_genes(fh):
    lines = []
    skip_tags = [ 'EC_number', 'gene' ]
    skip_features = [ 'REFERENCE', 'ncRNA' ]
    gene = []
    record = None
    for line in fh:
        line = line.rstrip()
        cols = line.split("\t")
        if line.startswith(">Feature"):
            if record and record["feature"] not in skip_features:
                gene.append(record)
            if len(gene) > 1: yield seqname, gene
            gene = []
            seqname = line.split()[1]
        elif len(cols) == 3:
            start, end, feature = cols
            if record and record["feature"] not in skip_features:
                gene.append(record)
            if feature == "gene":
                if len(gene) > 1: yield seqname, gene
                gene = []
            record = {
                "feature": feature,
                "desc": {},
                "coords": [ (start, end) ]
            }
        elif len(cols) == 2:
            start, end = cols
            record["coords"] += [ (start, end) ]
        elif len(cols) == 5:
            tag, value = cols[3:5]
            if tag not in skip_tags:
                if tag in record["desc"]:
                    if not isinstance(record["desc"][tag], list):
                        record["desc"][tag] = [ record["desc"][tag] ]
                    record["desc"][tag].append(value)
                else:
                    record["desc"][tag] = value
        else:
            raise ValueError('Unexpected line in .tbl: ' + line)
    if record:
        if record and record["feature"] not in skip_features:
            gene.append(record)
        if len(gene) > 1:
            yield seqname, gene

def fix_cds_coords(gene_record, cds_record):
    gene_start, gene_end = gene_record["coords"][0]
    start_partial = '<' in gene_start
    end_partial = '>' in gene_end
    if start_partial or end_partial:
        frame = int(cds_record["desc"]["codon_start"])
        g_start = int(gene_start.replace('<', ''))
        g_end   = int(gene_end.replace('>', ''))
        cds_start = cds_record["coords"][0][0]
        cds_end   = cds_record["coords"][-1][1]

        if start_partial:
            start = int(cds_start.replace('<', ''))
            diff = abs(start - g_start)
            if 0 < diff < 3:
                cds_record["coords"][0] = (gene_start, cds_record["coords"][0][1])
                cds_record["desc"]["codon_start"] = frame + diff
        if end_partial:
            end = int(cds_end.replace('>', ''))
            diff = abs(g_end - end)
            if 0 < diff < 3:
                cds_record["coords"][-1] = (cds_record["coords"][-1][0], gene_end)

def fix_short_intron(gene_record, child_record):
    if len(child_record["coords"]) > 1:
        start, end = child_record["coords"][0]
        start2, end2 = child_record["coords"][1]
        start = start.replace('<', '')
        exon_len = abs(int(end) - int(start)) + 1
        intron_len = abs(int(end) - int(start2)) - 1
        if exon_len <= 3:
            child_record["coords"].pop(0)
            child_record["coords"][0] = '<' + start2, end2
            gene_record["coords"][0] = '<' + start2, gene_record["coords"][0][1]
            if "codon_start" in child_record["desc"]:
                frame = int(child_record["desc"]["codon_start"])
                child_record["desc"]["codon_start"] = frame + 3 - exon_len

with open(tbl_file) as tbl:
    with open(tbl_out_file, 'w') as out:
        prev_seqname = ""
        for seqname, gene in get_genes(tbl):
            if seqname not in skip_sequences:
                lines = []
                if seqname != prev_seqname:
                    lines += [ ">Feature {seqname}".format(seqname = seqname) ]
                    lines += ncRNA[seqname]
                prev_seqname = seqname
                locus_tag = next(locus_tags)
                gene_record = gene[0]
                for record in gene:
                    if record["feature"] == "CDS":
                        fix_cds_coords(gene_record, record)
                    if record["feature"] == "CDS" or record["feature"] == "mRNA":
                        fix_short_intron(gene_record, record)
                for record in gene:                    
                    feature, desc, coords = record["feature"], record["desc"], record["coords"]
                    desc["locus_tag"] = locus_tag
                    lines += wrap_coords(feature, coords)
                    lines += wrap_desc(desc)
                for line in lines:
                    out.write(line + '\n')
