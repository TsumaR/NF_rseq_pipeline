#cd q/usr/bin/env python
#!/opt/conda/bin python

import sys
import argparse
import os
import pysam

from version import __version__

print(sys.path)

DEFAULT_OUT_BAM = 'out.bam'
TEMPORARY_SUFFIX = '.tmp'
TAG = 'gr'
TYPE = 'Z'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('bam', help='bam file, e.g. aln.bam')
    ap.add_argument('tag', help='region tag name, e.g. exon -> gr:Z:exon is added')
    ap.add_argument('-o', '--output-bam', dest='obam', action='store',
                    metavar='FILE', default=DEFAULT_OUT_BAM,
                    help='output bam (default=%s)' % DEFAULT_OUT_BAM)
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = ap.parse_args()

    fout = args.obam + TEMPORARY_SUFFIX if args.bam == args.obam else args.obam
    ibam = pysam.AlignmentFile(args.bam, 'rb')
    obam = pysam.AlignmentFile(fout, 'wb', template=ibam)
    for read in ibam:
        read.set_tag(TAG, args.tag, TYPE)
        obam.write(read)
    ibam.close()
    obam.close()

    if args.bam == args.obam:
        os.rename(fout, args.obam)

if __name__ == '__main__':
    main()
