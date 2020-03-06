import re

class Gene:
    def __init__(self, gtf):
        self.gene_id = gtf.attr.gene_id
        self.gene_name = gtf.attr.gene_name
        self.gene_biotype = gtf.attr.gene_biotype
        self.strand = gtf.strand
        self.chrom = gtf.seqname
        self.transcripts = []

    def add_transcript(self, tr):
        self.transcripts.append(tr)

class Transcript:
    def __init__(self, gtf):
        self.transcript_id = gtf.attr.transcript_id
        self.transcript_name = gtf.attr.transcript_name
        self.transcript_biotype = gtf.attr.transcript_biotype
        self.gene_id = gtf.attr.gene_id
        self.strand = gtf.strand
        self.chrom = gtf.seqname
        self.exon = []
        self.add_exon(gtf)

    def add_exon(self, gtf):
        self.exon.append(Exon(gtf))

    def coverage(self, bw):
        if bw.chroms(self.chrom) is None:
            []
        else:
            coordinate = []
            coverage = []
            for e in sorted(self.exon, key=lambda e: e.start):
                coordinate.extend(range(e.start, e.stop))
                coverage.extend(bw.values(self.chrom, e.start, e.stop))
            return (coordinate, coverage)

    def length(self):
        return sum(list(map(lambda e: e.length, self.exon)))

class Exon:
    def __init__(self, gtf):
        self.start = int(gtf.start) - 1
        self.stop = int(gtf.end)
        self.length = self.stop - self.start

class Record:
    HEADER = ('seqname', 'source', 'feature', 'start', 'end',
              'score', 'strand', 'frame')

    def __init__(self, line):
        for h in Record.HEADER:
            self.__dict__[h] = None
        self.attr = None
        self.__parse(line)

    def __parse(self, line):
        cols = line.split('\t')
        for i, h in enumerate(Record.HEADER):
            self.__dict__[h] = cols[i]
        self.attr = Attr(cols[8])

class Attr:
    ATTRS = ('gene_id', 'gene_name', 'gene_biotype', 'transcript_id',
             'transcript_name', 'transcript_biotype')

    def __init__(self, attr_str):
        for a in Attr.ATTRS:
            self.__dict__[a] = None
        self.__parse(attr_str)

    def __parse(self, attr_str):
        for pair in attr_str.split(';'):
            if pair.strip() != '':
                key, value = pair.strip().split(' ', 1)
                self.__dict__[key] = re.sub('"', '', value)

class Reader:
    @classmethod
    def each(self, fh):
        for line in fh:
            if not line.startswith('#'):
                yield Record(line.rstrip())
