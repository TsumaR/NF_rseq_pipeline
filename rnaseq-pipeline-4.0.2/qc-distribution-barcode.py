#!/usr/bin/env python

from collections import defaultdict
import argparse

import pysam

from lib import Extend, Qc
from version import __version__

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('bam', help='bam file')
    ap.add_argument('-o', '--output', dest='fout', action='store',
                    metavar='FILE', default=None,
                    help='output file, default=stdout')
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = ap.parse_args()

    hist_mis = defaultdict(Qc.BamMisInDelHist)
    hist_ins = defaultdict(Qc.BamMisInDelHist)
    hist_del = defaultdict(Qc.BamMisInDelHist)
    pos_mis = defaultdict(Qc.BamMisInDelPos)
    pos_ins = defaultdict(Qc.BamMisInDelPos)
    pos_del = defaultdict(Qc.BamMisInDelPos)
    with pysam.AlignmentFile(args.bam, 'r') as fh:
        for bam in fh:
            if bam.is_unmapped or bam.is_secondary:
                continue
            (query, bc, umi) = bam.query_name.split('_')
            cnt = Qc.BamMisInDel(bam)
            hist_mis[bc].add(cnt.n_mis)
            hist_ins[bc].add(cnt.n_ins)
            hist_del[bc].add(cnt.n_del)
            pos_mis[bc].add(cnt.pos_mis)
            pos_ins[bc].add(cnt.pos_ins)
            pos_del[bc].add(cnt.pos_del)

    with Extend.IO.sopen(args.fout) as fo:
        print('\t'.join(('barcode', 'data_type', 'x', 'y')), file=fo)
        for bc in hist_mis.keys():
            for v in zip(*(pos_mis[bc].x(), pos_mis[bc].y())):
                vv = tuple([bc] + list(v))
                print('%s\tpos_mismatch\t%d\t%d' % (vv), file=fo)
            for v in zip(*(pos_ins[bc].x(), pos_ins[bc].y())):
                vv = tuple([bc] + list(v))
                print('%s\tpos_insertion\t%d\t%d' % (vv), file=fo)
            for v in zip(*(pos_del[bc].x(), pos_del[bc].y())):
                vv = tuple([bc] + list(v))
                print('%s\tpos_deletion\t%d\t%d' % (vv), file=fo)
            for v in zip(*(hist_mis[bc].x(), hist_mis[bc].y())):
                vv = tuple([bc] + list(v))
                print('%s\tdist_mismatch\t%d\t%d' % (vv), file=fo)
            for v in zip(*(hist_ins[bc].x(), hist_ins[bc].y())):
                vv = tuple([bc] + list(v))
                print('%s\tdist_insertion\t%d\t%d' % (vv), file=fo)
            for v in zip(*(hist_del[bc].x(), hist_del[bc].y())):
                vv = tuple([bc] + list(v))
                print('%s\tdist_deletion\t%d\t%d' % (vv), file=fo)

if __name__ == '__main__':
    main()
