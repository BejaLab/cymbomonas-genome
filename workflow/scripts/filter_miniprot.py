
for input_file in snakefile.input:
    with open(input_file) as fh:
        for line in fh:
            if not line.startswith('#'):
                seqname, source, feature, start, end, score, strand, frame, attribute = line.split('\t')
                if feature == 'mRNA':
                    :
