
import re
from sys import stdin
from collections import defaultdict

attr_re = re.compile(r'([^=]+)=([^;]+);?')
upper_re = re.compile('[^A-Z]')
file_re = re.compile('file_\d+_|_')

def rename_id(ID, source):
    prefix = upper_re.sub('', source)[0:2]
    return prefix + file_re.sub('', ID)

tr_count = defaultdict(int)
for line in stdin:
    line = line.strip()
    if not line.startswith('#'):
        seqname, source, feature, start, end, score, strand, frame, attribute = line.split('\t')
        attrs = dict(attr_re.findall(attribute))
        if not line.endswith(';'):
            line += ';'
        if feature == 'transcript':
            attrs['Parent'] = attrs['locus'].replace('_', '')
            geneID = rename_id(attrs['geneID'], source)
            tr_count[geneID] += 1
            attrs['ID'] = '{locus}_{gene}.t{count}'.format(locus = attrs['Parent'], gene = geneID, count = tr_count[geneID])
            tr_id = attrs['ID']
        if feature == 'CDS' or feature == 'exon':
            attrs['Parent'] = tr_id
            attrs['ID'] = '{parent}.{feature}'.format(parent = attrs['Parent'], feature = feature)
        attribute = ''.join([ '{}={};'.format(key, value) for key, value in attrs.items()])
        line = '\t'.join([ seqname, source, feature, start, end, score, strand, frame, attribute ])
    print(line)
