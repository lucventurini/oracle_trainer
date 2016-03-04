#!/usr/bin/env python3
#coding: utf_8


import re
import sys,subprocess
from math import log
import io

class blastLine(object):
    '''This class provides the attributes for a typical blast (-m8) line. Attributes:
    query = the query id.
    subject = subject accession number/id
    identity = identity percentage.
    align_len = alignment length.
    mismatches = number of aligned bases with mismatch. (int)
    gaps = number of gaps.
    qStart = query Start.
    qEnd = query end.
    sStart = subject start.
    sEnd = subject end.
    e = evalue.
    bits: bit Score.'''
    def __init__(self,line):
        fields=line.rstrip().split('\t')
        self.query,self.subject=fields[0:2]
        self.identity=float(fields[2])
        self.align_len,self.mismatches,self.gaps,qStart,qEnd,sStart,sEnd=tuple(int(i) for i in fields[3:10])
        self.e=float(fields[10])
        self.bits=float(fields[11])
