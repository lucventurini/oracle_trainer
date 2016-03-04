#!/usr/bin/env python3


import sys,argparse
from scipy import mean

# Tue Dec  3 18:57:31 CET 2013

# **Summary Stats**
# Annotation:     /tmp/tmp4q45LN
# Predictions:                    /home/luca/Project_Melanzana/Annotation/Gene_annotation/NewMaker/final_alignment.maker.output/Training/Kcross/Test_0/predictions.dataset_0.gtf

# Gene Sensitivity                17.58%
# Gene Specificity                14.33%
# Transcript Sensitivity          18.01%
# Transcript Specificity          14.33%
# Exon Sensitivity                31.22%
# Exon Specificity                35.95%
# Nucleotide Sensitivity          91.09%
# Nucleotide Specificity          66.55%

def kcross(inps, f_score=False):

    '''This function parses the input eval files and calculates the mean spec., sens.
    at the gene, exon and nucleotide level.'''

    stats=dict()

    for inp in inps:
        record=False
        for line in open(inp):
            line=[x for x in line.rstrip().split() if x!='']
            if len(line)<3: continue
            if " ".join([line[0],line[1]])=="Gene Sensitivity":
                record=True
            if record:
                if " ".join(line[:2]) not in stats: stats[" ".join(line[:2])]=[]
                line[-1]=line[-1].rstrip("%")
                stats[" ".join(line[:2])].append(float(line[-1]))
            if " ".join([line[0],line[1]])=="Nucleotide Specificity": break

#print(stats, file=sys.stderr)

    for stat in sorted(stats.keys()):
        yield "\t".join(["{0}".format(stat), "{0}".format(mean(stats[stat]))])

def main():
    parser=argparse.ArgumentParser('''Quick utility to calculate the mean of multiple Eval files,
to obtain the final mean values for the cross-validation.''')
    parser.add_argument("eval", nargs='+', help="The list of eval files.")
    args=parser.parse_args()
    for line in kcross(args.eval): print(line)



if __name__=='__main__': main()
