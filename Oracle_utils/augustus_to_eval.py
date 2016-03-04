#!/usr/bin/env python3


import sys,argparse

def aug_to_eval(gtf):
    '''Main function. It receives in input a GTF file and yields the corrected lines (as a field *list*).'''
    
    if not isinstance(gtf, file):
        gtf=open(gtf)

    for line in gtf:
        line=line.rstrip()
        if  line=='' or line[0]=="#": continue
        if "a-posteriori" in line: break #END
        fields=line.split("\t")
        try:
            if fields[1]!="AUGUSTUS" or fields[2] in ("gene", "transcript"): continue
        except IndexError:
            continue
        fields[3]=int(fields[3])
        fields[4]=int(fields[4])
        seq,start_end=fields[0].split("_")
        start,end=start_end.split("-")
        fields[0]=seq
        start=int(start)
        fields[3]+=start-1
        fields[4]+=start-1
        yield fields

def main():
    
    parser=argparse.ArgumentParser("Script to convert Augustus GTF files into EVAL-compatible files.")
    parser.add_argument('gtf', type=argparse.FileType('r'))
    parser.add_argument('out', nargs='?', default=sys.stdout, type=argparse.FileType('w'))
    args=parser.parse_args()

    for line in aug_to_eval(args.gtf):
        print(*line, sep="\t", file=args.out)
    
if __name__=='__main__': main()
