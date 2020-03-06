from io import StringIO
import contextlib
import gzip
import os
import pandas as pd
import re
import sys
import zipfile

from . import Extend

class AnalysisEnv:
    def __init__(self):
        self.sh = '/bin/bash'
        self.ldir = 'log'
        self.sdir = 'summary'
        self.is_ajob = False
        self.is_bjob = False
        self.fqsub = None

class AbstractRNASeq(object):
    def __init__(self, smpl, conf, temp, env):
        self.smpl = smpl
        self.fconf = conf
        self.ftemp = temp
        self.env = env
        self.indent_lv = ' ' * 2
        if self.env.is_ajob:
            self.indent_lv = ''
            self.__output_func = self.__process_array_job
        else:
            self.__output_func = self.__process_for

    def __iprint(self, output, fo):
        print('{}{}'.format(self.indent_lv, output), file=fo)

    def __output_qsub(self, fo, fshell):
        with open(fshell) as fh:
            for line in fh:
                l = line.strip()
                if re.search('-pe', l):
                    l = '{} {}'.format(l, self.env.cpu)
                print(l, file=fo)

    def __output_shell(self, fo, fshell):
        with open(fshell) as fh:
            print(fh.read().strip(), file=fo)

    def __output_shell_with_indent(self, fo, fshell):
        with open(fshell) as fh:
            for line in fh:
                self.__iprint(line.strip(), fo)

    def __process_array_job(self, fo):
        print('cdir=$(pwd)', file=fo)
        print('i=$(( ${SGE_TASK_ID} - 1 ))', file=fo)
        print('id=${seqid_lst[$i]}', file=fo)
        print('cd $id', file=fo)
        self.__output_shell_with_indent(fo, self.ftemp)
        print('cd $cdir', file=fo)
        print('', file=fo)

    def __process_for(self, fo):
        print('cdir=$(pwd)', file=fo)
        print('for id in ${seqid_lst[@]}; do', file=fo)
        self.__iprint('cd $id', fo)
        self.__output_shell_with_indent(fo, self.ftemp)
        self.__iprint('cd $cdir', fo)
        print('done', file=fo)

    def _after_template(self, fo):
        print('#-- template end --#', file=fo)
        print('', file=fo)

    def _before_template(self, fo):
        print('#-- template start --#', file=fo)
        print('', file=fo)

    def build(self, fout):
        with Extend.IO.sopen(fout) as fo:
            #-- env --#
            print('#!{}'.format(self.env.sh), file=fo)
            if self.env.is_bjob or self.env.is_ajob:
                self.__output_qsub(fo, self.env.fqsub)
            if self.env.is_ajob:
                print('#$ -t 1:{}'.format(self.smpl.length), file=fo)
            print('', file=fo)
            #-- sequence id array --#
            print('seqid_lst=({})'.format(' '.join(self.smpl.df.sequence_id)), file=fo)
            print('', file=fo)
            #-- number of cores --#
            print('cpu={}'.format(self.env.cpu), file=fo)
            print('', file=fo)
            #-- config -#
            self.__output_shell(fo, self.fconf)
            #-- template --#
            self._before_template(fo)
            self.__output_func(fo)
            self._after_template(fo)
        Extend.OS.makedirs(self.env.ldir)

class Preparation:
    def __init__(self, smpl):
        self.smpl = Sample(smpl)

    def __mk_symlink(self, src, dst):
        if not os.path.exists(src):
            print('{} is not found.'.format(src), file=sys.stderr)
        elif os.path.exists(dst):
            print('{} is already existed.'.format(dst), file=sys.stderr)
        else:
            os.symlink(src, dst)

    def mk_symlinks(self, root=None, pref='raw'):
        if root is not None:
            root = os.path.abspath(root)
        cdir = os.getcwd()
        for i, item in self.smpl.df.iterrows():
            Extend.OS.makedirs(item.sequence_id)
            os.chdir(item.sequence_id)
            if self.smpl.is_paired:
                for i in (1, 2):
                    idx = 'fastq{}'.format(i)
                    dst = '{}_{}.fastq.gz'.format(pref, i)
                    src = item[idx] if root is None else '{}/{}'.format(root, item[idx])
                    self.__mk_symlink(src, dst)
            else:
                dst = '{}.fastq.gz'.format(pref)
                src = item.fastq1 if root is None else '{}/{}'.format(root, item.fastq1)
                self.__mk_symlink(src, dst)
            os.chdir(cdir)

class Sample:
    def __init__(self, smpl):
        self.df = pd.read_csv(smpl, sep='\t')
        self.is_paired = self.__check_paired()
        self.length = len(self.df.index)

    def __check_paired(self):
        if 'fastq1' in self.df.columns:
            return ('fastq2' in self.df.columns)
        else:
            raise pd.errors.ParserError('Invalid headers in sample file')

class SmartSeq(AbstractRNASeq):
    def __init__(self, smpl, conf, temp, env):
        super().__init__(smpl, conf, temp, env)

class Summary:
    def __init__(self, smpl, conf, temp, env):
        self.fsmpl = smpl
        self.fconf = conf
        self.ftemp = temp
        self.env = env

    def __output_shell(self, fo, fshell):
        with open(fshell) as fh:
            print(fh.read().strip(), file=fo)
        print('', file=fo)

    def __output_qsub(self, fo, fshell):
        with open(fshell) as fh:
            for line in fh:
                l = line.strip()
                if re.search('-pe', l):
                    l = '{} {}'.format(l, self.env.cpu)
                print(l, file=fo)

    def build(self, fout):
        with Extend.IO.sopen(fout) as fo:
            #-- env --#
            print('#!{}'.format(self.env.sh), file=fo)
            if self.env.is_bjob or self.env.is_ajob:
                self.__output_qsub(fo, self.env.fqsub)
            #-- number of cores --#
            print('cpu={}'.format(self.env.cpu), file=fo)
            print('', file=fo)
            #-- config -#
            self.__output_shell(fo, self.fconf)
            #-- summary dir, sample file --#
            print('summary_dir={}'.format(self.env.sdir), file=fo)
            print('smpl={}'.format(self.fsmpl), file=fo)
            print('', file=fo)
            #-- template --#
            self.__output_shell(fo, self.ftemp)
        Extend.OS.makedirs(self.env.sdir)

