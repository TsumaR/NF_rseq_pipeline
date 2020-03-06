#!/usr/bin/env python

import argparse
import math
import numpy as np
import pandas as pd

import pyBigWig

from lib import Extend, Gtf
from version import __version__

FISO = 'isoform_abund.tab'
LOWER = 10.0
BIN = 20

def coverage_bin(cov, seqlen, bins):
    bin_size = seqlen // bins + 1
    mod = bin_size * bins - seqlen
    start = 0
    stop = -1
    vals = []
    for i in range(0, bins):
        start = stop + 1
        stop = start + bin_size - (mod if i == 0 else 0) - 1
        c = list(filter(lambda v: not math.isnan(v), cov[1][start:(stop + 1)]))
        vals.append(np.mean(c) if len(c) > 0 else 0.0)
    return vals

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gtf', help='GTF file')
    ap.add_argument('bw', help='BigWig file')
    ap.add_argument('--isoform', dest='fiso', action='store',
                    metavar='FILE', default=FISO,
                    help='output file, default=' + FISO)
    ap.add_argument('--lower-limit', dest='low', action='store',
                    metavar='FLOAT', default=LOWER,
                    help='transcript lower limit, default=%.1f' % LOWER)
    ap.add_argument('--bin-size', dest='bins', action='store',
                    metavar='INT', default=BIN,
                    help='coverage into bin size, default=%d' % BIN)
    ap.add_argument('-o', '--output', dest='fout', action='store',
                    metavar='FILE', default=None,
                    help='output file, default=stdout')
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = ap.parse_args()

    tb = pd.read_csv(args.fiso, sep='\t')
    tar = tb[tb.TPM > args.low]
    tar_tr = {}
    for tid in tar.transcript_id:
        tar_tr[tid] = True

    trs = {}
    with open(args.gtf) as fh:
        for gtf in Gtf.Reader.each(fh):
            if gtf.attr.transcript_id not in tar_tr:
                continue
            if gtf.attr.transcript_id in trs:
                tr = trs[gtf.attr.transcript_id]
                tr.add_exon(gtf)
            else:
                tr = Gtf.Transcript(gtf)
                trs[tr.transcript_id] = tr

    with Extend.IO.sopen(args.fout) as fo:
        print('\t'.join(('transcript_id', 'strand', 'index', 'coverage')), file=fo)
        with pyBigWig.open(args.bw) as fh:
            for tr in trs.values():
                cov = tr.coverage(fh)
                vals = coverage_bin(cov, tr.length(), args.bins)
                if tr.strand == '-':
                    vals.reverse()
                for i, v in enumerate(vals):
                    print('%s\t%s\t%d\t%f' % (tr.transcript_id, tr.strand, i, v), file=fo)

if __name__ == '__main__':
    main()
