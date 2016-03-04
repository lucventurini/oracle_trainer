#!/usr/bin/env python3
#coding: utf_8


import re
import sys,subprocess
from math import log
import io

class Batch(object):
    '''This class functions as interface to parse the batch line produced by the batchBlast program.'''

    def __init__(self,handle): self.__handle=handle

    def __iter__(self): return self

    def __next__(self):
        bLine=batchLine()
        while True:
            line=self.__handle.readline()
            if line=='': raise StopIteration
            try: self.record=bline(line); break
            except:
                if line[0]=='#': return line
                bLine=batchLine()
        return self.record

class batchLine(object):
    '''This class defines the batch line. Arguments:
    query = the query id.
    subject = subject accession number
    identity = identity percentage.
    align_len = alignment length.
    e = e-value score.
    bits = bit score.
    qStart = query Start.
    qEnd = query end.
    sStart = subject start.
    sEnd = subject end.
    qLen = query length.
    sLen = subject length.
    qCov = query coverage.
    sCov = subject coverage.
    title = hit description.
    '''

    def __init__(self): pass

    def __call__(self,line):
        self.line=line
        fields=self.__line.rstrip().split('\t')
        #'Query','Subject','Identity','Al_len','E-val','BitScore','Query_start','Query_end','Sbj_start','Sbj_end','Subj_len','Query_cov','Subj_cov','GI_def'
        self.query=fields[0]
        self.subject=fields[1]
        self.identity=float(fields[2])
        self.align_len=int(fields[3])
        self.e=float(fields[4])
        self.bits=float(fields[5])
        self.qStart=int(fields[6])
        self.qEnd=int(fields[7])
        self.sStart=int(fields[8])
        self.sEnd=int(fields[9])
        self.qLen=int(fields[-5])
        self.sLen=int(fields[-4])
        self.qCov=float(fields[-3])
#        self.qLen=round(self.align_len*100/self.qCov,2)
        self.sCov=float(fields[-2])
        self.title=fields[-1]

    def __str__(self):
        try:
            fields=tuple(str(i) for i in (self.query,self.subject,self.identity,self.align_len,self.e,\
                self.bits,self.qStart,self.qEnd,self.sStart,self.qLen,self.sLen,self.qCov,self.sCov,self.title))
            return '\t'.join(fields)
        except NameError: raise NameError

