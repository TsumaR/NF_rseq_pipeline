#!/usr/bin/env python

import argparse
import os
import pandas as pd
import sys

from lib import RNASeq
from version import __version__

def get_analysis_env(args):
    env = RNASeq.AnalysisEnv()
    env.cpu = args.cpu
    env.is_ajob = args.ajob
    env.is_bjob = args.bjob
    if args.fqsub is None:
        env.fqsub = '{}/qsub_opt.sh'.format(os.path.abspath(os.path.dirname(__file__)))
    else:
        env.fqsub = args.fqsub
    return(env)

def cmd_prepare(args):
    prep = RNASeq.Preparation(args.sample)
    prep.mk_symlinks(args.root)

def cmd_smart(args):
    smpl = RNASeq.Sample(args.sample)
    ftemp = 't_smart_pe.sh' if smpl.is_paired else 't_smart_se.sh'
    ftemp = '{}/{}'.format(os.path.abspath(os.path.dirname(__file__)), ftemp)
    env = get_analysis_env(args)
    smart = RNASeq.SmartSeq(smpl, args.config, ftemp, env)
    smart.build(args.fout)

def cmd_summary(args):
    smpl = RNASeq.Sample(args.sample)
    ftemp = 't_summary_pe.sh' if smpl.is_paired else 't_summary_se.sh'
    ftemp = '{}/{}'.format(os.path.abspath(os.path.dirname(__file__)), ftemp)
    env = get_analysis_env(args)
    summ = RNASeq.Summary(args.sample, args.config, ftemp, env)
    summ.build(args.fout)

def common_opts(op):
    op.add_argument('-a', '--array-job', dest='ajob', action='store_true',
                    help='output shell script as array job')
    op.add_argument('-b', '--batch-job', dest='bjob', action='store_true',
                    help='output shell as batch job')
    op.add_argument('-o', '--output', dest='fout', action='store',
                    metavar='FILE', default=None,
                    help='output shell name, default=stdout')
    op.add_argument('-q', '--qsub-shell', dest='fqsub', action='store',
                    metavar='FILE', default=None,
                    help='qsub option shell, default=None')
    op.add_argument('-t', '--threads', dest='cpu', action='store',
                    metavar='INT', default=1, type=int,
                    help='number of threads, default=1')

def main():
    #-- global arguments --#
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    sp = ap.add_subparsers()

    #-- for preparation --#
    sp_p = sp.add_parser('prep', help='prepare to run rnaseq')
    sp_p.add_argument('sample', help='sample information file')
    sp_p.add_argument('-r', '--root-dir', dest='root', action='store',
                      metavar='DIR', default=None,
                      help='path to raw data root directory')
    sp_p.set_defaults(handler=cmd_prepare)

    #-- for smart-seq --#
    sp_s = sp.add_parser('smart', help='build shell script for SMART-Seq or Bead-seq')
    sp_s.add_argument('sample', help='sample information file')
    sp_s.add_argument('config', help='configuration shell file for SMART-Seq or Bead-seq')
    common_opts(sp_s)
    sp_s.set_defaults(handler=cmd_smart)

    #-- for summary --#
    sp_m = sp.add_parser('summary', help='summarize results')
    sp_m.add_argument('sample', help='sample information file')
    sp_m.add_argument('config', help='configuration shell file')
    common_opts(sp_m)
    sp_m.set_defaults(handler=cmd_summary)

    #-- parse arguments --#
    args = ap.parse_args()
    if hasattr(args, 'handler'):
        args.handler(args)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
