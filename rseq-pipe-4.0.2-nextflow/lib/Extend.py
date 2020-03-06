import contextlib
import gzip
import os
import sys

class IO:
    @classmethod
    @contextlib.contextmanager
    def sopen(self, fname=None):
        fh = open(fname, 'w') if fname and fname != '-' else sys.stdout
        try:
            yield fh
        finally:
            if fh is not sys.stdout:
                fh.close()

    @classmethod
    @contextlib.contextmanager
    def zopen(self, fname):
        fh = gzip.open(fname) if fname.endswith('.gz') else open(fname)
        try:
            yield fh
        finally:
            fh.close()

class OS:
    @classmethod
    def makedirs(self, path):
        if not os.path.exists(path):
            os.makedirs(path)
        else:
            print('%s is already existed.' % path, file=sys.stderr)

