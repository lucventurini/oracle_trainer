#!/usr/bin/env python3


import sys,argparse,gzip,re,tempfile,os
from operator import itemgetter, attrgetter
from Bio.Blast.NCBIXML import parse as xparser

from multiprocessing import Pool,Manager,Process
from multiprocessing import active_children as active
import time

'''This program will convert a Blast XML output into a conservation sequence.
The pseudocode for the procedure can be found in .
The algorithm functions as follows:
1- ConSeq = [unaligned]*len(sequence)
2- sort HSP by score
3- for i in xrange(1, length(sequence)):
   4- for each HSP in sorted(HSPs, order=score)
      5- if HSP extends to position i AND
      5b- ConSeq[i]==unaligned
          6- ConSeq[i]=HSP[i]
'''
def splitString(line, n=60):
    line=line.strip()
    length=len(line)
    stop_line=length-(length%60)
    for start in range(0,stop_line-1,60):
        yield line[start:start+60]
    yield line[stop_line:]


def conseq(record, symbol_table={}, from_file=False, transitions=False, gaps=False, unaligned=False, keep_symbols=False, with_name=True):
    '''Main function. It follows the pseudocode published by Korf et al. in 2001.
    The arguments regarding the type of ConSEQ are:
    - gaps: include gaps in the ConSEQ
    - transistions: include transistions in the ConSeq
    - unaligned: include unaligned positions in the ConSEQ
    - keep_symbols: use symbols instead of numbers for the ConSEQ
    Moreover, the following arguments are used to control the behavior of the function:
    - from_file: if set to True, it indicates that the "record" feature is an
                 XML mono-record Blast file which has to read internally.
    - symbol_table: table of the symbols, calculated by define_symbols.
                    If not defined, the relative function will be called internally.
    - with_name: If set to True (default) the output is a full FASTA file.
                 If set to False (e.g. in prepare_data) only '''

    if symbol_table=={}:
        symbol_table=define_symbols(unaligned=unaligned, transitions=transitions, gaps=gaps) #Calculate on the fly the symbol table

    if from_file:
        record_name=record[:]
        record=next(xparser(open(record)))

    transition={
        'A': {'G': 1, 'R': 1},
        'C': {'T': 1, 'Y': 1},
        'G': {'A': 1, 'R': 1},
        'T': {'C': 1, 'Y': 1},
        'R': {'A': 1, 'G': 1},
        'Y': {'C': 1, 'T': 1},
        'N': {},
        }

    sort_key="score" #Add configuration from args

    conseq_sequence=[None]*record.query_length
    my_descriptions=sorted(record.descriptions, key=attrgetter(sort_key))
    my_alignments=[]
    for description in my_descriptions:
        #Reorder alignments so that they match the order of descriptions
        my_alignments.append( next(filter(
                lambda alignment: alignment.accession==description.accession, record.alignments))) 

    for alignment in my_alignments:
        #Stop looking if all positions have already been defined
        if len([x for x in conseq_sequence if not x])==0: break
        for hsp in alignment.hsps:
            #Discard the information of gaps.
            correspondence=[corr for corr in zip( hsp.query,hsp.match,hsp.sbjct) if corr[0]!="-"]

            for position,corr in zip(range(hsp.query_start-1, hsp.query_end), correspondence):
#               print(position, corr, file=sys.stderr)
               #If the position had already been defined by another HSP, continue
               if conseq_sequence[position]!=None: continue 
               query_pos,match,sbjct_pos=corr
               token=None 
               if sbjct_pos=="-": token="-"
               elif match==" ":
                   if sbjct_pos in transition[query_pos]: token="/" #Is this a transition?
                   else: token=":" #If not, it's a mismatch
               else: token="|"
               conseq_sequence[position]=token
    
    final_conseq=[]
    for position in range(len(conseq_sequence)):
        if not conseq_sequence[position]:
            #Not defined positions are unaligned.
            conseq_sequence[position]="."
        if not transitions and conseq_sequence[position]=="/":
            conseq_sequence[position]=":" #Change transition to mismatch
        if not gaps and conseq_sequence[position]=="-":
            conseq_sequence[position]=":" #Change gap to mismatch
        if not unaligned and conseq_sequence[position]==".": 
            conseq_sequence[position]=":" #Change unaligned to mismatch
        token=conseq_sequence[position]
        if keep_symbols:
            final_conseq.append(token)
        else:
            final_conseq.append(symbol_table[token]) #Append the value

    assert len([x for x in final_conseq if x!=None and x!=''])==record.query_length
    #Format the list into a string
    final_conseq="".join([str(letter) for letter in final_conseq]) 
    lines=[line for line in splitString(final_conseq)]

    if from_file: os.remove(record_name)
    if with_name:
            lines="\n".join(lines)
            lines=">{0}\n{1}\n".format(record.query, lines)
    else:
        lines="".join(lines)

    return lines

def define_symbols(unaligned=False, gaps=False, transitions=False):
    '''This function formats the symbolic table used
    for the creation of the conseq file.'''
    symbols=[ ":", "|"]
    order={
        ':': 1, #Mismatch
        '|': 2, #Match
        '.': 3, #Unaligned
        '/': 4, #Transition
        '-': 5, #Gap
        }


    if unaligned: symbols.append(".")
    if gaps: symbols.append('-')
    if transitions: symbols.append("/")

    #Order the symbols in accordance with what we have to keep
    symbols=sorted(symbols, key=order.get)

    #Assign a progressive value to each retained symbol
    symbol_table=dict((symbol,i) for i,symbol in enumerate(symbols))
    return symbol_table

def wrapper(record, args, symbol_table={}, out=None, lock=None, from_file=False):
    '''This function reads the arguments from the argparse object
    and communicates them to the Conseq function. This structure makes it possible
    to call conseq from other programs (e.g. prepare_data)'''

    if out==None: out=sys.stdout

    conseq_sequence=conseq(record, transitions=args.transitions, 
                           gaps=args.gaps, unaligned=args.unaligned,
                           keep_symbols=args.keep_symbols,
                           symbol_table=symbol_table, from_file=from_file)
    lock.acquire()
    if out!=sys.stdout: out=open(out,'a')
    print(conseq_sequence, file=out, end="\n\n")
    if out!=sys.stdout: out.close()
    lock.release()
    return

def splitBlast(xml):
    header=''
    bool_header=True #Flag to understand whether I am still parsing the header
    temp='' #String with all the lines to print
    final_lines="  </BlastOutput_iterations>\n</BlastOutput>\n"
    tempXml=tempfile.NamedTemporaryFile(delete=False, suffix=".xml")
    for line in xml:
        if line.rstrip().lstrip()=="<Iteration>":
            #Finished parsing one seq results, time to print
            if bool_header:
                bool_header=False #Finished parsing header
                temp+=line
            else:
                print(header, temp, final_lines, end='', file=tempXml)
                temp=''
                tempXml.close()
                yield tempXml.name
                temp+=line

            tempXml=tempfile.NamedTemporaryFile(delete=False, suffix=".xml")

        else:
            if bool_header: header+=line
            else: temp+=line

    print(header,temp, end='', file=tempXml) #Do NOT include final_lines in the last output
    tempXml.close()
    yield tempXml.name
            

def main():
    parser=argparse.ArgumentParser("This utility calculates a ConSeq file as defined by the Brent Lab.")
    parser.add_argument('-x', "--keep_symbols", default=False, action="store_true",
                        help="Keep the symbols of the original ConSeq, not their numerical representation.")
    parser.add_argument("-u", "--unaligned", default=False, action="store_true",
                        help="Keep also information on unaligned portions.")
    parser.add_argument('--processors', default=1, type=int, help="Number of processors to use. Default: 1")
    parser.add_argument('--debug', default=False, action='store_true', help="Flag. If set, multiprocessing will be disabled to allow finer debugging.")
    parser.add_argument("-t", "--transitions", default=False, action="store_true",
                        help="Keep also information on transitions.")
    parser.add_argument("-g", "--gaps", default=False, action="store_true",
                        help="Keep also information on alignment gaps.")
    parser.add_argument('--out', default=sys.stdout, type=argparse.FileType('w'),
                        help="Output file. Default: stdout")
    parser.add_argument('xml', type=argparse.FileType('r'), help="The XML Blast file to analyze.")
    args=parser.parse_args()

    if re.search("\.gz$", args.xml.name): args.xml=gzip.open(args.xml.name) #If the file is gzip-compressed, open it as such
    if args.out!=sys.stdout:
        args.out.close()
        args.out=args.out.name

    pool=Pool(processes=args.processors)
 #The principal process
    manager=Manager()
    lock=manager.RLock()
    symbols=define_symbols(unaligned=args.unaligned, gaps=args.gaps, transitions=args.transitions)
    #principal=active()[0]
    # parser=xparser(args.xml)
    jobs=[]
    if args.out==sys.stdout: args.out=None
    for temporary_file in splitBlast(args.xml):
        # while True:
        #     if len(active())-1<args.processors: break #Rest while no processors are available. active() will make any finished processes join.
        if args.debug:
            wrapper(temporary_file, args, symbol_table=symbols, out=args.out, lock=lock, from_file=True)
        else:
                   
            job=pool.apply_async(wrapper, args=(temporary_file, args),
                                 kwds={'symbol_table':symbols,
                                       'out': args.out,
                                       'lock': lock,
                                       'from_file': True})
            jobs.append(job)

    if not args.debug:
        for job in jobs: job.get()
    pool.close()
    pool.join()
                             

    return


if __name__=='__main__': main()
