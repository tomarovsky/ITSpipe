#!/usr/bin/env python3
import re
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
    infile = open(args.input, 'r').readlines()
    outfile = open(args.output, 'a')
    pattern_len = len(args.pattern)
    outfile.write("".join(infile[:4])) # header
    reverse = False # flag
    prev_tlen = None
    for line in infile[4:]:
        line = line.strip().split("\t")
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line[:11]
        bitwise_flags = '\t'.join(line[11:])
        if seq.startswith(args.pattern) and pos == args.reference_start and int(tlen) >= 0 and reverse is False:
            tlen = str(int(tlen) - pattern_len)
            seq = seq[pattern_len:]
            qual = qual[pattern_len:]
            cigar = cigar_left_trimmer(cigar, pattern_len)
            pos = str(int(pos) + pattern_len)
            reverse = True
            prev_tlen = int(pos)
        elif reverse and int(tlen) == 0 - prev_tlen:
            tlen = str(int(tlen) + pattern_len)
            seq = seq[pattern_len:]
            qual = qual[pattern_len:]
            cigar = cigar_left_trimmer(cigar, pattern_len)
            pos = str(int(pos) - pattern_len)
            reverse = False
        line = "\t".join([qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, bitwise_flags])
        outfile.write(line)
        outfile.write("\n")
    outfile.close()


if __name__ == '__main__':
    parser = ArgumentParser(description="Script for trimming the SAM file from the high coverage area")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help="input SAM file")
    group_required.add_argument('-o', '--output', type=str, help="output trimmed SAM file")
    group_required.add_argument('-p', '--pattern', type=str, help="sequence with high coverage")
    group_required.add_argument('-s', '--reference_start', type=str, help="1-based leftmost coordinate")
    args = parser.parse_args()
    main()