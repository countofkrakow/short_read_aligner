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
    with open(fname, 'r') as f:
        read_id = f.readline()
        while read_id:
            read_sequence = f.readline()
            ret.append(re.sub(r'[^AGTC]', 'T', read_sequence.strip()))

            # ignore this
            f.readline()
            f.readline()
            read_id = f.readline()
    return ret


