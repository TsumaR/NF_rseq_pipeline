#!/usr/bin/env python

from collections import defaultdict
import argparse

import pysam

from lib import Extend
from version import __version__

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('bam', help='bam file')
    ap.add_argument('-o', '--output', dest='fout', action='store',
                    metavar='FILE', default=None,
                    help='output file, default=stdout')
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = ap.parse_args()

    isize = defaultdict(int)
    with pysam.AlignmentFile(args.bam, 'r') as fh:
        for bam in fh:
            if bam.is_unmapped or bam.is_secondary:
                continue
            if bam.is_proper_pair and bam.is_read1:
                tlen = bam.tlen if bam.tlen > 0 else bam.tlen * -1
                isize[tlen] += 1

    with Extend.IO.sopen(args.fout) as fo:
        print('\t'.join(('insert_size', 'count')), file=fo)
        for s in sorted(isize.keys()):
            print('%d\t%d' % (s, isize[s]), file=fo)

if __name__ == '__main__':
    main()
