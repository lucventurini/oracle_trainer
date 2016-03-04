#!/usr/bin/env python3


import sys
import re
import os
import subprocess
import tempfile
import time
from myRecords import GFF
from Bio import SeqIO, Seq, SeqRecord
from .markov import Markov
from numpy import arange
from operator import itemgetter
import io

##In the maker-derived GTF the stop codon is EXTERNAL to any other annotation, be it CDS or UTR.
##Start codons are instead considered part of the CDS.

########################################

def get_eval(evalFile):
    '''This function evaluates an EVAL output file and outputs it F1 scores (gene, exon, nucleotide)'''
    record=False
    stats={}
    if not isinstance(evalFile, io.IOBase): evalFile=open(evalFile)

    for line in evalFile:
        line=[x for x in line.rstrip().split() if x!='']
        if len(line)<3: continue
        if " ".join([line[0],line[1]])=="Gene Sensitivity": record=True
        if record:
            # if " ".join(line[:2]) not in stats: stats[" ".join(line[:2])]=[]
            line[-1]=line[-1].rstrip("%")
            stats[" ".join(line[:2])]=float(line[-1])
        if " ".join([line[0],line[1]])=="Nucleotide Specificity": break

    gene_stats=[stats[key]/100 for key in [key for key in list(stats.keys()) if "gene" in key.lower()]]
    exon_stats=[stats[key]/100 for key in [key for key in list(stats.keys()) if "exon" in key.lower()]]
    nucl_stats=[stats[key]/100 for key in [key for key in list(stats.keys()) if "nucl" in key.lower()]]

    # if sum(gene_stats)>0: gene_f1=2*( prod(gene_stats)/sum(gene_stats))
    # else: gene_f1=0

    # if sum(exon_stats)>0: exon_f1=2*( prod(exon_stats)/sum(exon_stats))
    # else: exon_f1=0

    # if sum(nucl_stats)>0: nucl_f1=2*(prod(nucl_stats)/sum(nucl_stats))
    # else: nucl_f1=0

    return gene_stats + exon_stats + nucl_stats

########################################

def parseval(parseval_file):

    '''Function to parse very simply the parseval output and recover the F1 scores.'''

    if isinstance(parseval, io.IOBase): par=parseval_file
    else: par=open(parseval_file)

    gene=None
    exon=None
    nucl=None
    vals={'gene': 0,
          'exon': 0,
          'nucl': 0}
         

    while True:
        try: line=next(par)
        except StopIteration: break
        if line[0]=="=" or line.rstrip().lstrip()=="": val=None #Avoid blank or comment lines, and re-instantiate the dictionary
        elif "CDS structure comparison" in line: val='gene'
        elif "Exon structure comparison" in line: val='exon'
        elif "Nucleotide-level comparison" in line: val='nucl'
        elif "F1 Score" in line:
            fields=[x for x in line.rstrip().split() if x not in ('',"F1")]
            
            if len(fields)==1: fields=[x for x in re.split("\.\.+", line.rstrip()) if x!=''] #Split eliminating those pesky multiple dots
            try: vals[val]=float(fields[1])
            except ValueError:
                if fields[1]=="--": vals[val]=0
                else: raise ValueError("could not convert string to float: {0}".format(fields))

    return vals['gene'], vals['exon'], vals['nucl']

######################################

def slimGff(gff):

    '''This function removes spurious tags from the GeneID GFF3 files.
    It yields the corrected lines one by one.'''


    #Check that the input is not a file/buffer
    if isinstance(gff, io.IOBase) or isinstance(gff, io.BufferedReader):
        gff=gff
    else: gff=GFF.GFF3(open(gff))

    for line in gff:
        if isinstance(line, bytes): line=line.decode('UTF-8') #Necessary to decode the lines from bytestring to str
        if line[0]=="#": yield line

        else:
            line=GFF.gffLine(line)
            attrs={}
            for att in line.attributes:
                if att not in ("ID","Parent", "Name"): continue
                attrs[att]=line.attributes[att]
            line.attributes=attrs
            yield str(line)

#####################################

def fasta_dict(fasta):

    '''Dictionary to store the information regarding the number of times we see a sequence'''

    dictionary={}
    for seq in SeqIO.parse(open(fasta),'fasta'):
        seq=str(seq.seq)
        dictionary[seq]=dictionary.get(seq,0)+1
    return dictionary

######################################

def gff3_to_gtf(gff, gtf, sequence, lock=None):
    '''This function takes as input two filehandles and a SeqIO.index:
    - gff (read), the input
    - gtf (write), the output
    - seqIndex, the SeqIO.index object used to access the sequences.
    It will convert the input into a fixed gtf file using the gff3_to_eval_gtf and validate_gtf utilities.'''

    gtf=open(gtf, 'w')

    seqIndex=SeqIO.index(sequence, 'fasta')

    conversion=subprocess.Popen([os.path.join(os.environ['EVAL_GTF'], 'gff3_to_gtf.pl'), gff.name], shell=False, stdout=subprocess.PIPE,
                                stderr=open(os.devnull, 'w') )

    currSeq=None
    tempGtf=tempfile.NamedTemporaryFile(suffix=".gtf", mode="w+") #STUPID DAMNED ****ING EVALUATE. It truncates if it sees a gtf. It is bound to happen sometimes :-(
    prefix=re.sub("\.gtf$", "", tempGtf.name) #rstrip breaks things. Bad rstrip.

    for line in conversion.stdout:
        line=line.decode('UTF-8').rstrip()
        if line=='': continue
        chrom=line.split()[0]
        if chrom!=currSeq:
#            print("In the cycle", file=sys.stderr)
            if currSeq:
                tempFa=tempfile.NamedTemporaryFile(mode="w+")
                print(seqIndex[currSeq].format('fasta'), file=tempFa)
                tempFa.flush()
                tempGtf.flush()
                subprocess.call(['/opt/eval-2.2.8/validate_gtf.pl', '-f',
                                 tempGtf.name, tempFa.name], shell=False,
                                stdout=open(os.devnull,'w'), stderr=open(os.devnull,'w')) #I do not want any stupid output clogging the logs
                tempFa.close()
                for line in open(prefix+".fixed.gtf"): print(line, file=gtf)
                tempGtf.close()
                os.remove(prefix+".fixed.gtf")
                tempGtf=tempfile.NamedTemporaryFile(suffix=".gtf")
                prefix=re.sub("\.gtf$", "", tempGtf.name)
            currSeq=chrom
        print(line, file=tempGtf)

    if lock: lock.acquire()
    tempFa=tempfile.NamedTemporaryFile()
    if lock: lock.release()
    try:
        print(seqIndex[currSeq].format('fasta'), file=tempFa)
    except: pass
    tempFa.flush()
    tempGtf.flush()
    subprocess.call([os.path.join(os.environ['EVAL_GTF'],'validate_gtf.pl'), '-f',
                     tempGtf.name, tempFa.name], shell=False,
                    stdout=open(os.devnull,'w'), stderr=open(os.devnull,'w'))
    tempFa.close()
    for line in open(prefix+".fixed.gtf"):
        line=line.rstrip()
        if line=="": continue #Avoid blank lines
        print(line, file=gtf)
    tempGtf.close()
    os.remove(prefix+".fixed.gtf")
    gff.close()
    tempFa.close()
    tempGtf.close()
    gtf.close()
    return gtf

####################################

def doGeneID(gtf, param, seqfile, lock=None):
    '''This function will create a fixed GTF from the GFF3 GID prediction. This will allow for quick comparison using EVAL.'''

#    seqIndex=SeqIO.index(seqfile.name, 'fasta')
    if lock: lock.acquire()
    gff=tempfile.NamedTemporaryFile(mode="w+")
    if lock: lock.release()
#    print(gff.name, file=sys.stderr)
    geneid=subprocess.Popen(["geneid", "-3P", param, seqfile], shell=False, stdout=subprocess.PIPE)
    lineCount=0

    first=True

    for line in slimGff(geneid.stdout): #we need to slim down the gff
        print(line, file=gff)
    gff.flush()
    gtf=gff3_to_gtf(gff, gtf, seqfile)
    gtf.close()
    if len(open(gtf.name).__next__().split("\t"))==9: return True
    else: return False
    

###################################


def optimize(records, args, paramFile, lock=None, pool=None):
    '''This function checks the various OWF and EWF possible factors and finds the best one. Inspiration: Optimization_new.sh'''

    dataFolder=os.path.dirname(paramFile.name)
    scores=[]
    if args.debug: bDel=False
    else: bDel=True

    optimFolder=os.path.join(dataFolder, "GID_Optimization")
    if not os.path.exists(optimFolder): os.makedirs(optimFolder)

    report=open("/".join([dataFolder,'optimization_report.txt']),'w')
    header=['sizeFactor', 'exonWeight', 'Sn_gene', 'Sp_gene','Sn_exon', 'Sp_exon', 'Sn_nucl','Sp_nucl']
    print(*header, file=report)
    for sizeFactor in arange(args.iOwf, args.fOwf+args.dOwf, args.dOwf): #We have to use arange as (x)range does NOT support floating points.
        for exonWeight in arange(args.iEwf, args.fEwf+args.dEwf, args.dEwf): 
            param=os.path.join(optimFolder, "geneid_sF_{0}_eF_{1}.param".format(sizeFactor,exonWeight))
            if not os.path.exists(param) or os.stat(param).st_size==0:
                param=open(param,'w')
                for line in writeParam(paramFile, sizeFactor, exonWeight, size_exon_diff=args.size_exon_diff): print(line, file=param)
                param.close()
                param=param.name

            refList=open(
                os.path.join(optimFolder, "geneid_sF_{0}_eF_{1}.ref.list".format(sizeFactor,exonWeight)),
                'w')

            predList=open(
                os.path.join(optimFolder, "geneid_sF_{0}_eF_{1}.pred.list".format(sizeFactor,exonWeight)),
                'w')
            
            gtf_folder=os.path.join(optimFolder, "geneid_sF_{0}_eF_{1}".format(sizeFactor,exonWeight))
            if not os.path.exists(gtf_folder): os.makedirs(gtf_folder)
            param_jobs=[]
            for record in records:
                name, folder=record 
                gtf=os.path.join(gtf_folder, "{0}.gtf".format(name))
                if (not os.path.exists(gtf)) or (os.stat(gtf).st_size==0):
                    fasta=os.path.join(folder,"{0}.masked.fa".format(name))
                    assert os.path.exists(fasta)
                    if pool!=None:
                        param_jobs.append(pool.apply_async(doGeneID, args=(gtf, param, fasta)))
                    else:
                        doGeneID(gtf, param, fasta) #Has GID found any models for the sequence?
                print(os.path.abspath(gtf), file=predList)
                print(os.path.join(folder, "{0}.noutr.gtf".format(name)), file=refList)

                # else:
                #     print(sizeFactor, exonWeight, file=sys.stderr)
                #     if len(open(gtf).next().split("\t"))==9: hits.append(name) #Has GID found any models for the sequence?

            for job in param_jobs: job.get()

            evalFile=os.path.join(optimFolder, "geneid_sF_{0}_eF_{1}.eval".format(sizeFactor, exonWeight))
            refList.close()
            predList.close()              

            if (not os.path.exists(evalFile)) or (os.stat(evalFile).st_size==0):
                evalFile=open(evalFile,'w')
                ##This somewhat complicated writing out is necessary as evaluate_gtf.pl uses the *order* of the files to guess the matching.
                evaluation=subprocess.call([os.path.join(os.environ['EVAL_GTF'],'evaluate_gtf.pl'), '-Aq',
                                            refList.name, predList.name], stdout=evalFile, stderr=open(os.devnull,'w'))
            else:
                evalFile=open(evalFile)

            gene_sn, gene_sp, exon_sn, exon_sp, nucl_sn, nucl_sp = get_eval(evalFile.name)
#            print(gene_f1, exon_f1, nucl_f1, file=sys.stderr)
            
            score={'sizeFactor': sizeFactor,
                   'exonWeight': exonWeight,
                   'Sn_gene': round(gene_sn,5),
                   'Sp_gene': round(gene_sp,5),
                   'Sn_exon': round(exon_sn,5),
                   'Sp_exon': round(exon_sp,5),
                   'Sn_nucl': round(nucl_sn,5),
                   'Sp_nucl': round(nucl_sp,5)
                   }
            scores.append(score)
            line=[]
            for key in header:
                line.append(score[key])
#            print(*line, sep="\t", file=report)
#            param.close()
    
    bestVals=sorted(scores, key=itemgetter('Sn_gene', 'Sp_gene', 'Sn_exon', 'Sp_exon', 'Sn_nucl', 'Sp_nucl'), reverse=True)
    for val in bestVals:
        line=[val[key] for key in header]
        print(*line, file=report, sep="\t")

    chosen=bestVals[0]

    for line in writeParam(paramFile, chosen['sizeFactor'], chosen['exonWeight'], args.size_exon_diff):
        yield line #The function in the main program will write everything out

##############################################

def prepare_files(records, folder, args, keep=True):

    '''This function recovers the necessary files from the folders
    created by the prepare_data program.'''

#    out=open(out.name,'w')


    starts=open(os.path.join(folder, "starts.fasta"),'w')
    false_starts=open(os.path.join(folder, "false_starts.fasta"),'w')

    donors=open(os.path.join(folder, "donors.fasta"),'w')
    false_donors=open(os.path.join(folder, "false_donors.fasta"),'w')

    acceptors=open(os.path.join(folder, "acceptors.fasta"),'w')
    false_acceptors=open(os.path.join(folder, "false_acceptors.fasta"),'w')

    cds=open(os.path.join(folder, "cds.fasta"),'w')
    introns=open(os.path.join(folder, "introns.fasta"),'w')

    startCount=0
    donorCount=0
    acceptorCount=0

    if args.masked: addition=".masked"
    else: addition=''
    master_folder=folder[:]

    for record,folder in records:

        #Recover CDS
        try:
            print(SeqIO.read(open(
                        os.path.join(folder,"{0}.cds{1}.fa".format(record, addition))),'fasta').format('fasta'),
                  file=cds, end='')
        except: raise ValueError(os.path.join(folder,"{0}.cds{1}.fa".format(record, addition)))

        #Recover Start
        print(SeqIO.read(open(
                    os.path.join(folder,"{0}.start{1}.fa".format(record, addition))),'fasta').format('fasta'),
              file=starts, end='')

        #Recover introns
        for intron in SeqIO.parse(open(
                os.path.join(folder,"{0}.introns{1}.fa".format(record, addition) )),'fasta'):
            print(intron.format('fasta'), file=introns, end='')

        #Recover donors
        for donor in SeqIO.parse(open(
                os.path.join(folder,"{0}.donors{1}.fa".format(record, addition) )),'fasta'):
            print(donor.format('fasta'), file=donors, end='')
        

        #Recover Acceptors                        
        for acceptor in SeqIO.parse(open(
                os.path.join(folder,"{0}.acceptors{1}.fa".format(record, addition))),'fasta'):
            print(acceptor.format('fasta'), file=acceptors, end='')

        #Recover False Starts
        for line in open(os.path.join(folder,"{0}.false_starts{1}.fa".format(record, addition) )):
            startCount+=1
            line=line.rstrip()
            seq=SeqRecord.SeqRecord(Seq.Seq(line), id=str(startCount), description='')
            print(seq.format('fasta'), file=false_starts, end='')
                        
        for line in open(os.path.join(folder,"{0}.false_donors{1}.fa".format(record, addition) )):
            donorCount+=1
            line=line.rstrip()
            seq=SeqRecord.SeqRecord(Seq.Seq(line), id=str(donorCount), description='')
            print(seq.format('fasta'), file=false_donors, end='')

        for line in open(os.path.join(folder,"{0}.false_acceptors{1}.fa".format(record, addition))):
            acceptorCount+=1
            line=line.rstrip()
            seq=SeqRecord.SeqRecord(Seq.Seq(line), id=str(acceptorCount), description='')
            print(seq.format('fasta'), file=false_acceptors, end='')

    print(time.ctime(), "Finished retrieving files for {0}".format(master_folder))


def prepare_template(records, folder, args, keep=True):

    '''This function prepares the unoptimized.param file
    by calling markov.Markov on all the data sets (cds 
    vs. introns, starts vs. false_starts, etc.). If the
    neecessary files are not present it calls prepare_files.'''

    #Start Signal

    starts=os.path.join(folder, "starts.fasta")
    false_starts=os.path.join(folder, "false_starts.fasta")
    start_freqs=os.path.join(folder, "start_freqs.txt")
    print(start_freqs, file=sys.stderr)
    start_markov=None

    if os.path.exists(start_freqs):
       try:
           start_markov=Markov(signal_probs = start_freqs, signal="Start", folder=folder,
                            signal_maximum=args.start_maximum, min_zscore=args.min_zscore,
                            cutoff=args.start_cutoff,
                            window_size = args.window_size, retrieve=True)
       except ValueError as error:
           if error.message=="Pre-calculated signals are greater than the requested maximum!":
               os.remove(start_freqs)
               start_reindex=True
           else:
               raise

    if start_markov==None:
        try:
            startsFasta=fasta_dict(starts)
        except:
            prepare_files(records, folder, args)
            try: startsFasta=fasta_dict(starts)
            except: raise
        false_startsFasta=fasta_dict(false_starts)
        start_markov=Markov(observed=startsFasta, background=false_startsFasta, 
                            signal="Start", folder=folder,
                            signal_maximum=args.start_maximum, min_zscore=args.min_zscore, 
                            cutoff=args.start_cutoff,
                            window_size = args.window_size)
        del startsFasta
        del false_startsFasta

    if not keep:
        os.remove(starts.name)
        os.remove(false_starts.name)
    
    #Donor Signal
    donors=os.path.join(folder, "donors.fasta")
    false_donors=os.path.join(folder, "false_donors.fasta")
    donor_freqs=os.path.join(folder, "donor_freqs.txt")
    donor_markov=None
    if os.path.exists(donor_freqs):
        try:
            donor_markov = Markov(signal_probs = donor_freqs, signal="Donor", folder=folder,
                                  signal_maximum=args.signal_maximum, min_zscore=args.min_zscore,
                                  cutoff=args.donor_cutoff,
                                  window_size = args.window_size, retrieve=True)
        except ValueError as error:
            if error.message=="Pre-calculated signals are greater than the requested maximum!":
                os.remove(donor_freqs)
            else:
                raise
    if donor_markov==None:
        donorsFasta=fasta_dict(donors)
        false_donorsFasta=fasta_dict(false_donors)
        donor_markov=Markov(observed=donorsFasta, background=false_donorsFasta, 
                            signal="Donor", folder=folder,
                            signal_maximum=args.signal_maximum, min_zscore=args.min_zscore, 
                            cutoff=args.donor_cutoff,
                            window_size = args.window_size)
        del donorsFasta
        del false_donorsFasta

    if not keep:
        os.remove(donors)
        os.remove(false_donors)
    
    #Acceptor signal
    acceptors=os.path.join(folder, "acceptors.fasta")
    false_acceptors=os.path.join(folder, "false_acceptors.fasta")
    acceptor_freqs = os.path.join(folder, "acceptor_freqs.txt")
    acceptor_markov=None
    if os.path.exists(acceptor_freqs):
        try:
            acceptor_markov = Markov(signal_probs = acceptor_freqs, signal="Acceptor", folder=folder,
                                     signal_maximum=args.signal_maximum, min_zscore=args.min_zscore,
                                     cutoff=args.acceptor_cutoff,
                                     window_size = args.window_size, retrieve=True)
        except ValueError as error:
            if error.message=="Pre-calculated signals are greater than the requested maximum!":
                os.remove(acceptor_freqs)
                acceptor_reindex=True
            else:
                raise

    if acceptor_markov==None:
        acceptorsFasta=fasta_dict(acceptors)
        false_acceptorsFasta=fasta_dict(false_acceptors)
        acceptor_markov=Markov(observed=acceptorsFasta, background=false_acceptorsFasta, 
                               signal="Acceptor", folder=folder,
                               signal_maximum=args.signal_maximum, min_zscore=args.min_zscore, 
                               cutoff=args.acceptor_cutoff,
                               window_size = args.window_size)
        del acceptorsFasta
        del false_acceptorsFasta

    if not keep:
        os.remove(acceptors)
        os.remove(false_acceptors)

    #Coding signals
    cds=os.path.join(folder, "cds.fasta")
    introns=os.path.join(folder, "introns.fasta")
    cds_initial = os.path.join(folder, "cds_initial_freqs.txt")
    cds_transition = os.path.join(folder, "cds_transition_freqs.txt")
    #Do not recalculate everything if it is already present
    cds_markov=None
    if os.path.exists(cds_initial) and os.path.exists(cds_transition):
        try:
            cds_markov = Markov( initial=cds_initial, transition=cds_transition, cds_maximum = args.cds_maximum, folder = folder, retrieve=True)
        except ValueError as error:
            if error.message=="Pre-calculated CDS kmers are greater than the requested maximum!": #we have increased the maximum information we want .. gotta recalculate
                os.remove(cds_initial)
                os.remove(cds_transition)
            else:
                raise
    if cds_markov==None:
        cdsFasta=fasta_dict(cds)
        intronsFasta=fasta_dict(introns)
        cds_markov=Markov(observed=cdsFasta, background=intronsFasta, cds_maximum=args.cds_maximum, folder=folder)
        
        del intronsFasta
        del cdsFasta

    if not keep:
        os.remove(cds)
        os.remove(introns)

    for line in open(args.template):
        if line.rstrip()=="Start_profile":
            for MMline in start_markov.print_markov():
                yield MMline
                #print(MMline, file=out)
        elif line.rstrip()=="Donor_profile":
            yield ""
            #print("", file=out)
            for MMline in donor_markov.print_markov():
                yield MMline
                #print(MMline, file=out)
        elif line.rstrip()=="Acceptor_profile":
            yield ""
            #print("", file=out)
            for MMline in acceptor_markov.print_markov(): 
                yield MMline
                #print(MMline, file=out)

        elif line.rstrip()=="Markov_oligo_logs_file":
            yield line.rstrip()
            #print(line.rstrip(), file=out)
            for MMline in cds_markov.print_markov():
                yield MMline
                #print(MMline, file=out)
        else:
            yield line.rstrip()
            #print(line.rstrip(), file=out)

def writeParam(original, sizeFactor, exonWeight, size_exon_diff=0):

    '''This function has the purpose to insert the sizeFactor and exonWeight at their right place in the geneid_param file.
    Arguments:
    original = the original parameter file
    sizeFactor = the SZ parameter
    exonWeight = the EW parameter
    size_exon_diff = the difference between sizeFactor and exonFactor. Default: 0.
    In tomato and potato, this value is 0.2'''

    lineCount=0
    original=open(original.name)
    for line in original:
        lineCount+=1
        if lineCount==27:
#            yield " ".join([str(i) for i in [sizeFactor+0.2]*4])
             yield " ".join([str(i) for i in [sizeFactor+size_exon_diff]*4])
        elif lineCount==30:
            yield " ".join([str(i) for i in [sizeFactor]*4])
        elif lineCount==36:
            yield " ".join([str(i) for i in [exonWeight]*4])
        else:
            yield line.rstrip()
