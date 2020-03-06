from collections import defaultdict
import numpy as np
import os
import re
import zipfile

class FastqcBasicStatistics:
    def __init__(self, rec, key):
        self.title = key
        self.gc_content = float(rec['%GC'])
        self.seqlen_range = rec['Sequence length']
        self.file_type = rec['File type']
        self.total_seq = int(rec['Total Sequences'])
        self.encoding = rec['Encoding']
        self.file_name = rec['Filename']
        seqlen = list(map(lambda v: int(v), self.seqlen_range.split('-')))
        self.seqlen_max = max(seqlen)

class FastqcBasicTable:
    def __init__(self, rec, key, x, y):
        self.title = key
        self.xname = x
        self.yname = y
        self.x = list(map(lambda v: np.mean(list(map(lambda vv: int(vv), v.split('-')))), rec[x]))
        self.y = list(map(lambda v: float(v), rec[y]))

class FastqcStats:
    def __init__(self, fzip):
        self.fzip = fzip
        self.basic_info = None
        self.base_qual = None
        self.qual = None
        self.gc = None
        self.seqlen = None
        self.__set(fzip)

    def __set(self, fzip):
        for key, rec in self.__each_data(fzip):
            if key == 'Basic Statistics':
                self.basic_info = FastqcBasicStatistics(rec, key)
            elif key == 'Per base sequence quality':
                self.base_qual = FastqcBasicTable(rec, key, 'Base', 'Mean')
            elif key == 'Per sequence quality scores':
                self.qual = FastqcBasicTable(rec, key, 'Quality', 'Count')
            elif key == 'Per sequence GC content':
                self.gc = FastqcBasicTable(rec, key, 'GC Content', 'Count')
            elif key == 'Sequence Length Distribution':
                self.seqlen = FastqcBasicTable(rec, key, 'Length', 'Count')

    def __each_data(self, fzip, fdat='fastqc_data.txt'):
        with zipfile.ZipFile(fzip, 'r') as zf:
            flst = zf.namelist()
            for f in filter(lambda f: os.path.basename(f) == fdat, flst):
                with zf.open(f) as fh:
                    key = None; data = []; header = []
                    for line in fh:
                        line = line.decode('utf-8').rstrip()
                        m = re.match('>>(.+?)\t', line)
                        if m is not None:
                            key = m.group(1)
                            continue
                        m = re.match('#', line)
                        if m is not None:
                            header = line.replace('#', '').split('\t')
                            continue
                        m = re.match('>>END_MODULE', line)
                        if m is not None:
                            rec = {}
                            if key == 'Basic Statistics':
                                for v in data:
                                    rec[v[0]] = v[1]
                            else:
                                tdata = list(zip(*data))
                                for i, h in enumerate(header):
                                    rec[h] = tdata[i]
                            yield(key, rec)
                            key = None; data = []; header = []
                            continue
                        data.append(line.split('\t'))

class BamMisInDel:
    def __init__(self, bam):
        self.n_mis = 0
        self.n_ins = 0
        self.n_del = 0
        self.pos_mis = {}
        self.pos_ins = {}
        self.pos_del = {}
        self.__parse(bam)

    def __parse(self, bam):
        aln = bam.get_aligned_pairs(with_seq=True)
        qidx = 0
        for i in range(bam.qstart, bam.qend):
            v = aln[i]
            if v[2] is not None and v[2].islower():
                self.pos_mis[v[0]] = 1
            elif v[1] is None:
                self.pos_ins[v[0]] = 1
            elif v[0] is None and v[2] is not None:
                self.pos_del[qidx] = 1
            if v[0] is not None:
                qidx = v[0]
        self.n_mis = len(self.pos_mis.keys())
        self.n_ins = len(self.pos_ins.keys())
        self.n_del = len(self.pos_del.keys())

class BamMisInDelHist:
    def __init__(self):
        self.data = defaultdict(int)

    def add(self, val):
        self.data[val] += 1

    def x(self):
        return sorted(self.data.keys())

    def y(self):
        return list(map(lambda k: self.data[k], self.x()))

class BamMisInDelPos:
    def __init__(self):
        self.data = defaultdict(int)

    def add(self, d):
        for k, v in d.items():
            self.data[k] += 1

    def x(self):
        return sorted(self.data.keys())

    def y(self):
        return list(map(lambda k: self.data[k], self.x()))

