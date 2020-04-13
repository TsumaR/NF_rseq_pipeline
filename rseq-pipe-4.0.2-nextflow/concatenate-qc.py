#!/usr/bin/env python

import argparse
import pandas as pd
import os
import sys

from version import __version__


def mk_table(qc, sid, read=None):
    tb = pd.read_csv(qc, sep='\t')
    if read is not None:
        tb['read'] = read
    tbex = tb.assign(Sample_ID = sid)
    colname = ['Sample_ID']
    colname.extend(tb.columns.values)
    return tbex[colname]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('smpl', help='sample information file')
    ap.add_argument('out_dir', help='output directory')
    ap.add_argument('-pe', dest='fpe', action='store', help='type of input file', default=None)
    ap.add_argument('-o', '--output', dest='fout', action='store',
                    metavar='FILE', default=None,
                    help='output file path, default=stdout')
    ap.add_argument('-t', '--type', dest='ptype', action='store')
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = ap.parse_args()

    # read input
    smpl = pd.read_csv(args.smpl, sep='\t')
    p_type = args.ptype
    LR = args.fpe

    data = None
    for i,sid in enumerate(smpl.Sample_ID):

        tb = None
        if LR != "LR":
            tb = mk_table(os.path.join(args.out_dir,  f'{smpl.iat[i,1]}', '08_qc', f'{sid}_{p_type}.txt'), sid)
        else:
            tb1 = mk_table(os.path.join(args.out_dir, f'{smpl.iat[i,1]}', '08_qc', f'{sid}_R1_{p_type}.txt'), sid, 'read1')
            tb2 = mk_table(os.path.join(args.out_dir, f'{smpl.iat[i,1]}', '08_qc', f'{sid}_R2_{p_type}.txt'), sid, 'read2')
            tb = pd.concat([tb1, tb2])
            
        if data is None:
            data = tb
        else:
            data = pd.concat([data, tb])

    if args.fout is not None:
        data.to_csv(args.fout, sep='\t', index=False)
    else:
        data.to_csv(sys.stdout, sep='\t', index=False)

if __name__ == '__main__':
    main()
