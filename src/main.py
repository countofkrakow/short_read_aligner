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
            return suffix_array[suffix_arr_idx]
    if l > r:
        return suffix_array[l]
    else:
        return suffix_array[l - 1]

def process_reads(reads, data, suffix_array, mismatches_allowed, min_mismatches_to_display):
    # process first hiv reads
    table = {}
    longest_prefix = ''
    longest_prefix_start = None
    accepted_alignments = 0
    for read_seq in reads:
        idx = binary_search(data, read_seq, suffix_array)
        prefix, mismatches = compare_read(data, idx, read_seq, mismatches_allowed)

        # check to see if there was an acceptable number of mismatches
        if mismatches:
            accepted_alignments += 1
            start, length = prefix
            if length > len(longest_prefix):
                longest_prefix = data[start:start+length]
                longest_prefix_start = start
            for mismatch in mismatches:
                if mismatch not in table:
                    table[mismatch] = 0
                table[mismatch] += 1


    print('Number of accepted alignments: %d' % accepted_alignments)
    print('Longest prefix: ' + longest_prefix)
    print('Start: ' + str(longest_prefix_start))

    sum_mismatches = 0
    for mismatches in table:
        count = table[mismatches]
        pos, nuc = mismatches
        if count >= min_mismatches_to_display:
            print('%d \t%s\t%d' % (pos, nuc, count))
        sum_mismatches += count

    print("Total number of mismatches: %d" % sum_mismatches)

    return table

def compare_read(data, read_idx, read_seq, max_mismatches):
    mismatches = []
    prefix = (read_idx, len(read_seq))
    for i, read_nucleotide in enumerate(read_seq):
        if len(mismatches) > max_mismatches:
            return prefix, []

        data_nucleotide = data[read_idx + i]
        if read_nucleotide != data_nucleotide:
            if not mismatches:
                # start position in genome and length
                prefix = (read_idx, i)
            mismatches.append((read_idx + i, read_nucleotide))
    return prefix, mismatches

if __name__ == '__main__':
    MAX_MISMATCHES_ALLOWED = 10
    READ_SEQUENCE_LENGTH = 100

    data = parse_fasta('hiv1.fasta')
    v1 = parse_fastq('visit1-hiv.fastq')
    v6 = parse_fastq('visit6-hiv.fastq')

    suffix_array = list(range(len(data) - READ_SEQUENCE_LENGTH))
    suffix_array.sort(key=lambda x: data[x:x + READ_SEQUENCE_LENGTH])

    print("Visit 1 of the HIV genome")
    v1_table = process_reads(v1, data, suffix_array, MAX_MISMATCHES_ALLOWED, 5)
    print()
    print("Visit 6 of the HIV genome")
    v6_table = process_reads(v6, data, suffix_array, MAX_MISMATCHES_ALLOWED, 5)
    print()

    v1_read_positions = set([a[0] for a in v1_table])
    v6_read_positions = set([a[0] for a in v6_table])

    v1_and_v6 = v1_read_positions & v6_read_positions
    in_v1_not_v6 = v1_read_positions - v6_read_positions
    print("Read positions in both v1 and v6 (%d)" % len(v1_and_v6))
    for pos in v1_and_v6:
        print('\t%d' % pos)

    print("Read positions in v1 and not in v6 (%d)" % len(in_v1_not_v6))
    for pos in in_v1_not_v6:
        print('\t%d' % pos)


