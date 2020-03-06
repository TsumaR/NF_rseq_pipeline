#!/usr/bin/env python

from collections import defaultdict
import argparse

import pysam

from lib import Extend, Qc
from version import __version__

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('qc_raw', help='raw fastqc zip')
    ap.add_argument('bc2', help='barcode fastq')
    ap.add_argument('aln', help='alignment bam')
    ap.add_argument('-o', '--output', dest='fout', action='store',
                    metavar='FILE', default=None,
                    help='output file, default=stdout')
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = ap.parse_args()

    qc_raw = Qc.FastqcStats(args.qc_raw)

    barcode = defaultdict(int)
    with Extend.IO.zopen(args.bc2) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                e1 = line.rstrip().split(' ')
                e2 = e1[0].split('_')
                barcode[e2[1]] += 1

    bcsum = sum(barcode.values())
    total = qc_raw.basic_info.total_seq
    unknown = total - bcsum

    region = defaultdict(lambda: defaultdict(int))
    with pysam.AlignmentFile(args.aln, 'r') as fh:
        for bam in fh:
            if bam.is_unmapped or bam.is_secondary:
                continue
            (query_name, bc, umi) = bam.query_name.split('_')
            region[bc][bam.get_tag('gr')] += 1

    with Extend.IO.sopen(args.fout) as fo:
        print('\t'.join(('barcode', 'alignment', 'count')), file=fo)
        print('unknown\tremoved\t%d' % unknown, file=fo)
        for bc, cnt in sorted(barcode.items(), key=lambda x: x[0]):
            mapped = reduce(lambda x, y: x + y, region[bc].values())
            unmapped = cnt - mapped
            print('%s\t%s\t%d' % (bc, 'unmapped', unmapped), file=fo)
            for k, v in region[bc].items():
                print('%s\t%s\t%d' % (bc, k, v), file=fo)

if __name__ == '__main__':
    main()
