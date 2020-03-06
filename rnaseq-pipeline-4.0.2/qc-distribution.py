#!/usr/bin/env python

from collections import defaultdict
import argparse

import pysam

from lib import Extend, Qc
from version import __version__

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cln', help='cleaned fastq')
    ap.add_argument('bam', help='bam file')
    ap.add_argument('--read2', dest='r2', action='store_true', default=False,
                    help='calculate about read2 of paired-end')
    ap.add_argument('-o', '--output', dest='fout', action='store',
                    metavar='FILE', default=None,
                    help='output file, default=stdout')
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = ap.parse_args()

    qc_cln = Qc.FastqcStats(args.cln)

    hist_mis = Qc.BamMisInDelHist()
    hist_ins = Qc.BamMisInDelHist()
    hist_del = Qc.BamMisInDelHist()
    pos_mis = Qc.BamMisInDelPos()
    pos_ins = Qc.BamMisInDelPos()
    pos_del = Qc.BamMisInDelPos()
    with pysam.AlignmentFile(args.bam, 'r') as fh:
        for bam in fh:
            if bam.is_unmapped or bam.is_secondary:
                continue
            if args.r2 and bam.is_read1:
                continue
            if not args.r2 and bam.is_read2:
                continue
            cnt = Qc.BamMisInDel(bam)
            hist_mis.add(cnt.n_mis)
            hist_ins.add(cnt.n_ins)
            hist_del.add(cnt.n_del)
            pos_mis.add(cnt.pos_mis)
            pos_ins.add(cnt.pos_ins)
            pos_del.add(cnt.pos_del)

    with Extend.IO.sopen(args.fout) as fo:
        print('\t'.join(('data_type', 'x', 'y')), file=fo)
        for v in zip(*(qc_cln.base_qual.x, qc_cln.base_qual.y)):
            print('pos_score\t%d\t%f' % v, file=fo)
        for v in zip(*(pos_mis.x(), pos_mis.y())):
            print('pos_mismatch\t%d\t%d' % v, file=fo)
        for v in zip(*(pos_ins.x(), pos_ins.y())):
            print('pos_insertion\t%d\t%d' % v, file=fo)
        for v in zip(*(pos_del.x(), pos_del.y())):
            print('pos_deletion\t%d\t%d' % v, file=fo)
        for v in zip(*(qc_cln.qual.x, qc_cln.qual.y)):
            print('dist_score\t%d\t%d' % v, file=fo)
        for v in zip(*(qc_cln.gc.x, qc_cln.gc.y)):
            print('dist_gc\t%d\t%d' % v, file=fo)
        for v in zip(*(qc_cln.seqlen.x, qc_cln.seqlen.y)):
            print('dist_seqlen\t%.1f\t%d' % v, file=fo)
        for v in zip(*(hist_mis.x(), hist_mis.y())):
            print('dist_mismatch\t%d\t%d' % v, file=fo)
        for v in zip(*(hist_ins.x(), hist_ins.y())):
            print('dist_insertion\t%d\t%d' % v, file=fo)
        for v in zip(*(hist_del.x(), hist_del.y())):
            print('dist_deletion\t%d\t%d' % v, file=fo)

if __name__ == '__main__':
    main()
