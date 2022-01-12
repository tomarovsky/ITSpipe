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
    outfile.write("".join(infile[:4])) # write header
    infile_without_header = infile[4:]

    for i in range(0, len(infile_without_header), 2):
        forward = infile_without_header[i].strip().split("\t")
        reverse = infile_without_header[i + 1].strip().split("\t")
        f_qname, f_flag, f_rname, f_pos, f_mapq, f_cigar, f_rnext, f_pnext, f_tlen, f_seq, f_qual = forward[:11]
        f_bitwise_flags = '\t'.join(forward[11:])
        r_qname, r_flag, r_rname, r_pos, r_mapq, r_cigar, r_rnext, r_pnext, r_tlen, r_seq, r_qual = reverse[:11]
        r_bitwise_flags = '\t'.join(reverse[11:])

        # if f_pos == args.pos and f_seq.startswith(args.pattern) and int(f_tlen) > 0:
        #     f_seq = f_seq[pattern_len:]
        #     f_qual = f_qual[pattern_len:]
        #     f_cigar = cigar_left_trimmer(f_cigar, pattern_len)
        #     r_seq = r_seq[pattern_len:]
        #     r_qual = r_qual[pattern_len:]
        #     r_cigar = cigar_left_trimmer(r_cigar, pattern_len)
        #     f_tlen = str(int(f_tlen) - pattern_len)  # +
        #     r_tlen = str(int(r_tlen) + pattern_len)
        #     f_pos = str(int(f_pos) + pattern_len)
        #     r_pos = str(int(r_pos) - pattern_len)
        #     r_pnext = f_pos
        #     f_pnext = r_pos
        if r_pos == args.pos and r_seq.startswith(args.pattern) and int(r_tlen) > 0:
            f_seq = f_seq[pattern_len:]
            f_qual = f_qual[pattern_len:]
            f_cigar = cigar_left_trimmer(f_cigar, pattern_len)
            r_seq = r_seq[pattern_len:]
            r_qual = r_qual[pattern_len:]
            r_cigar = cigar_left_trimmer(r_cigar, pattern_len)
            f_tlen = str(int(f_tlen) + pattern_len)
            r_tlen = str(int(r_tlen) - pattern_len)
            f_pos = str(int(f_pos) + pattern_len) #
            r_pos = str(int(r_pos) - pattern_len) #
            r_pnext = f_pos
            f_pnext = r_pos

        forward = "\t".join([f_qname, f_flag, f_rname, f_pos, f_mapq, f_cigar, f_rnext, f_pnext, f_tlen, f_seq, f_qual, f_bitwise_flags])
        reverse = "\t".join([r_qname, r_flag, r_rname, r_pos, r_mapq, r_cigar, r_rnext, r_pnext, r_tlen, r_seq, r_qual, r_bitwise_flags])
        outfile.write("%s\n" % forward)
        outfile.write("%s\n" % reverse)
    outfile.close()


if __name__ == '__main__':
    parser = ArgumentParser(description="Script for trimming the SAM file from the high coverage area")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help="input SAM file")
    group_required.add_argument('-o', '--output', type=str, help="output trimmed SAM file")
    group_required.add_argument('--pattern', type=str, help="sequence with high coverage")
    group_required.add_argument('--pos', type=str, help="1-based leftmost coordinate")
    args = parser.parse_args()
    main()