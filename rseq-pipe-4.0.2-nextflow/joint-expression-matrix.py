#!/usr/bin/env python

#version 4.0.3
#for Nextflow pipeline

import argparse
import pandas as pd
import sys

from lib import Extend, Gtf
from version import __version__

INFILE = 'gene_abund.tab'
ROOT_DIR = '.'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gtf', help='GTF file')
    ap.add_argument('smpl', help='sample information file')
    ap.add_argument('-i', '--input', dest='fin', action='store',
                    metavar='DIR', default=None,
                    help='gene abund file directory name, default=%s' % INFILE)
    ap.add_argument('--root', dest='root', action='store',
                    metavar='DIR', default=ROOT_DIR,
                    help='analysis root directory, default=%s' % ROOT_DIR)
    ap.add_argument('-o', '--output', dest='fout', action='store',
                    metavar='FILE', default=None,
                    help='output file path, default=stdout')
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = ap.parse_args()

    smpl = pd.read_csv(args.smpl, sep='\t')
    gene = {}
    with open(args.gtf) as fh:
        for gtf in Gtf.Reader.each(fh):
            if gtf.attr.gene_id not in gene:
                gene[gtf.attr.gene_id] = Gtf.Gene(gtf)
    genes = sorted(gene.values(), key=lambda v: v.gene_id)

    data = []
    for s in smpl.Sample_ID:
        print(s, file=sys.stderr)
        spath = '%s/%s' % (args.fin, f"{s}_R1_abund.tab")
        dat = pd.read_csv(spath, index_col='Gene ID', sep='\t')
        grp = dat.groupby(level=0)
        summary = grp.sum()
        data.append(summary.to_dict(orient='index'))

    with Extend.IO.sopen(args.fout) as fo:
        header = ['gene_id', 'chrom', 'gene_name', 'gene_biotype']
        header.extend(smpl.Sample_ID)
        print('\t'.join(header), file=fo)
        for rec in sorted(gene.values(), key=lambda v: v.gene_id):
            g = rec.gene_id
            vals= [rec.gene_id, rec.chrom, rec.gene_name, rec.gene_biotype]
            vals.extend(map(lambda v: v[g]['TPM'] if g in v else 0, data))
            print('\t'.join(map(lambda v: str(v), vals)), file=fo)

if __name__ == '__main__':
    main()
