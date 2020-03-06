#!/usr/bin/env python

from collections import defaultdict
import argparse

import pysam

from lib import Extend, Qc
from version import __version__

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('qc_raw', help='raw fastqc zip')
    ap.add_argument('qc_cln', help='cleaned fastqc zip')
    ap.add_argument('bam', help='bam file')
    ap.add_argument('--read2', dest='r2', action='store_true', default=False,
                    help='calculate about read2 of paired-end')
    ap.add_argument('-o', '--output', dest='fout', action='store',
                    metavar='FILE', default=None,
                    help='output file, default=stdout')
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = ap.parse_args()

    qc_raw = Qc.FastqcStats(args.qc_raw)
    qc_cln = Qc.FastqcStats(args.qc_cln)

    region = defaultdict(int)
    with pysam.AlignmentFile(args.bam, 'r') as fh:
        for bam in fh:
            if bam.is_unmapped or bam.is_secondary:
                continue
            if args.r2 and bam.is_read1:
                continue
            if not args.r2 and bam.is_read2:
                continue
            region[bam.get_tag('gr')] += 1

    mapped = sum(region.values())
    removed = qc_raw.basic_info.total_seq - qc_cln.basic_info.total_seq
    unmapped = qc_cln.basic_info.total_seq - mapped

    with Extend.IO.sopen(args.fout) as fo:
        print('\t'.join(('alignment', 'count')), file=fo)
        print('removed\t%d' % removed, file=fo)
        print('unmapped\t%d' % unmapped, file=fo)
        for k in ('intergenic', 'intron', 'ercc', 'exon'):
            print('%s\t%d' % (k, region[k] if k in region else 0), file=fo)

if __name__ == '__main__':
    main()
