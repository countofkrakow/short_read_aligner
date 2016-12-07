import re

def parse_fasta(fname):
    gene = ''
    for line in open(fname, 'r'):
        if '>' in line:
            if gene:
                return gene
        else:
            gene += re.sub(r'[^AGTC]', 'T', line.strip())
    return gene

def parse_fastq(fname):
    ret = []
    f = open(fname, 'r')
