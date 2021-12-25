#!/usr/bin/env python3
import re
import pysam
from argparse import ArgumentParser


def cigar_left_trimmer(cigar_line, pattern_len):
    cigar_matches = [{'type':m[1], 'length':int(m[0])} for m in re.findall(r'(\d+)([A-Z]{1})', cigar_line)]
    symbols = [l['type'] for l in cigar_matches]
    lengths = [l['length'] for l in cigar_matches]
    trimmed_cigar = ''
    for i in range(len(lengths)):
        if sum(lengths[:i+1]) > pattern_len:
            lengths[i] = sum(lengths[:i+1]) - pattern_len
            lengths, symbols = lengths[i:], symbols[i:]
            for l, s in zip(lengths, symbols):
                trimmed_cigar += ''.join([str(l), str(s)])
            break
        else:
            continue
    return trimmed_cigar

def main():
    infile = pysam.AlignmentFile(args.input, 'rb')
    outfile = pysam.AlignmentFile(args.output, 'wb', template=infile)
    pattern_len = len(args.pattern)

    for line in infile:
        if line.query_sequence.startswith(args.pattern) and line.reference_start == args.reference_start:
            line.template_length = line.template_length - pattern_len
            q = line.query_qualities
            line.query_sequence = line.query_sequence[pattern_len:]    # from the beginning of the sequence
            line.query_qualities = q[pattern_len:]
            line.cigarstring = cigar_left_trimmer(line.cigarstring, pattern_len)
            line.reference_start = line.reference_start + pattern_len
        outfile.write(line)
    outfile.close()

if __name__ == '__main__':
    parser = ArgumentParser(description="Script for trimming the BAM file from the high coverage area")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help="input BAM file")
    group_required.add_argument('-o', '--output', type=str, help="output BAM file")
    group_required.add_argument('-p', '--pattern', type=str, help="sequence with high coverage")
    group_required.add_argument('-s', '--reference_start', type=int, help="0-based leftmost coordinate")
    args = parser.parse_args()
    main()