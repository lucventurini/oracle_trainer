#!/usr/bin/env python3


import sys
from operator import attrgetter
from myRecords.GFF import GFF3

class gffRecord(object):
    #It might be ingenious to codify it as an object rather than using procedural programming..
    def __init__(self):
        self.sequences=[]

    def add(self, line):
        if line.feature=="CDS": self.sequences.append(line)       

    def printRecord(self):
        self.sequences=sorted(self.sequences, key=attrgetter('start'))
        for num in range(len(self.sequences)):
            feature=None
            sequence=self.sequences[num]
            if len(self.sequences)==1: feature="Single"
            else:
                if num==0:
                    if sequence.strand=="+"  or not sequence.strand: feature="Initial"
                    elif sequence.strand=="-": feature="Terminal"
                elif num==len(self.sequences)-1:
                    if sequence.strand=="+" or not sequence.strand: feature="Terminal"
                    elif sequence.strand=="-": feature="Initial"
                else:
                    feature="Internal"

            line=[sequence.chrom,
                  sequence.source,
                  feature,
                  sequence.start,
                  sequence.end,
                  sequence.score if sequence.score else ".",
                  sequence.strand if sequence.strand else ".",
                  sequence.phase if sequence.phase else ".",
                  sequence.parent]
            line=[str(el) for el in line]
            yield "\t".join(line)

def convert(gff):
    current=gffRecord()

    for line in GFF3(gff):
        if line.feature=="gene":
            if current.sequences: 
                for line in current.printRecord(): yield line
            current=gffRecord()
        else: current.add(line)

    if current.sequences:
        for line in current.printRecord(): yield line


if __name__=="__main__":

    for line in convert(sys.argv[1]): print(line)
