#!/usr/bin/env python3


import sys,re

def zff2gtf(zff):

    '''Script to convert ZFF format (SNAP) into a normal GTF. Source: http://www.bioperl.org/wiki/ZFF'''
    if isinstance(zff, file): inp=zff
    else:
        assert isinstance(zff, str)
        inp=open(zff)

    currSeq=None
    currModel=None
    for line in inp:
        stopLine=None
        if line[0]==">":
            currSeq=line.rstrip().lstrip(">")
            yield ""
        else:
            fields=line.rstrip().split("\t")
            feature,start,stop,strand,score,overhang5,overhang3,frame,model=fields
            start,stop=int(start),int(stop)
            score=float(score)
            if currModel!=model: currModel=model
            if feature in ("Esngl", "Einit"):
                if strand=="-":
                    initStart=stop-2
                    initStop=stop
                else:
                    initStart=start
                    initStop=start+2

                startLine=[currSeq, "snap", "start_codon",
                           initStart, initStop, 0, strand,
                           0, "gene_id \"{0}\"; transcript_id \"{0}\";".format(model)]

                startLine=[str(el) for el in startLine]
                yield "\t".join(startLine)
            if feature in ("Esngl", "Eterm"):
                if strand=="-":
                    termStart=start
                    termStop=start+2
                    start=start+3 #We have to move the stop codon outside the CDS
                else:
                    termStart=stop-2
                    termStop=stop
                    stop=stop-3
                stopLine=[currSeq, "snap", "stop_codon",
                           termStart, termStop, 0, strand,
                           overhang5, "gene_id \"{0}\"; transcript_id \"{0}\";".format(model)]
                stopLine=[str(el) for el in stopLine]

            line=[currSeq, "snap", "CDS",
                  start,stop,score,strand,
                  overhang5, "gene_id \"{0}\"; transcript_id \"{0}\";".format(model)]
            line=[str(el) for el in line]
            yield "\t".join(line)
            if stopLine: yield "\t".join(stopLine)

if __name__=='__main__':

    for line in zff2gtf(sys.argv[1]): print(line)
