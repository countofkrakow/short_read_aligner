from data import parse_fasta, parse_fastq


def create_index(data, substring_len=4):
    if len(data) < substring_len:
        raise IndexError('Substring length (%d) can\'t be longer than the data (%d).' % (substring_len, len(data)))

    idx = {}
    for i in range(len(data) - substring_len):
        substring = data[i:i+substring_len]
        if substring not in idx:
            idx[substring] = []
        idx[substring].append(i)
    return idx

def binary_search(data, read_seq, suffix_array):
    l, r = 0, len(suffix_array)
    while l < r:
        suffix_arr_idx = l + (r - l) // 2
        data_idx = suffix_array[suffix_arr_idx]
        seq_in_data = data[data_idx:data_idx+len(read_seq)]
        if seq_in_data < read_seq:
            l = suffix_arr_idx + 1
        elif seq_in_data > read_seq:
            r = suffix_arr_idx - 1
        else:
            return suffix_arr_idx
    if l > r:
        return l
    else:
        return l - 1

if __name__ == '__main__':
    MAX_MISMATCHES_ALLOWED = 10
    READ_SEQUENCE_LENGTH = 100

    hiv_genome = parse_fasta('hiv1.fasta')
    v1 = parse_fastq('visit1-hiv.fastq')
    v6 = parse_fastq('visit6-hiv.fastq')

    suffix_array = range(len(hiv_genome) - READ_SEQUENCE_LENGTH)
    suffix_array.sort(key=lambda x: hiv_genome[x:x+READ_SEQUENCE_LENGTH])


