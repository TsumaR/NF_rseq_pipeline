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
    tbex = tb.assign(sequence_id = sid)
    colname = ['sequence_id']
    colname.extend(tb.columns.values)
    return tbex[colname]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('smpl', help='sample information file')
    ap.add_argument('qc1', help='qc1 file')
    ap.add_argument('qc2', help='qc2 file for paired-end', nargs='?')
    ap.add_argument('-o', '--output', dest='fout', action='store',
                    metavar='FILE', default=None,
                    help='output file path, default=stdout')
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = ap.parse_args()

    # read input
    smpl = pd.read_csv(args.smpl, sep='\t')

    data = None
    for sid in smpl.sequence_id:
        tb = None
        if args.qc2 is None:
            tb = mk_table(os.path.join(sid, args.qc1), sid)
        else:
            tb1 = mk_table(os.path.join(sid, args.qc1), sid, 'read1')
            tb2 = mk_table(os.path.join(sid, args.qc2), sid, 'read2')
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
