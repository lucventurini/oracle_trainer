#!/usr/bin/env python3

import sys
import os
import re
import subprocess
import tempfile
import argparse
from multiprocessing import Pool, Manager, cpu_count
from Bio import SeqIO
from copy import copy
from Oracle_utils.kcross_mean import kcross
import textwrap
from Oracle_utils import geneid


def twinScan(folder, record, zhmm, outFolder, args, lock, add_conseq=False, add_estseq=False, ref_gtf=None):

        '''This function performs twinscan on the *masked* sequence for a given record.
        According to the values provided with the keywords
             - add_conseq
             - add_estseq
        twinscan will use (or not) the conseq and estseq files.'''
        lock.acquire()
#        print(folder, record, zhmm, outFolder, file=sys.stderr)
        lock.release()

        twinscan=args.twinscan
        alias=args.test_alias

        fasta=os.path.join(folder, record+".masked.fa")
        commandline=[twinscan]
        # if add_estseq:
        #       commandline.append("-pe")
        commandline+=[ zhmm, fasta ]
        out=os.path.abspath(os.path.join(outFolder, "{0}.gtf".format(record)))
        if os.path.exists(out) and os.stat(out).st_size:
                if args.force:
                        os.remove(out)
                else:
                        return out, 0

        out=open(out,'w')

        if add_conseq:
                if args.nscan:
                        align_name=os.path.join(folder, "{0}.{1}.align".format(record,alias))
                        commandline.append("-a={0}".format(align_name))

                else:
                        conseq=tempfile.NamedTemporaryFile(mode='w')
                        first=True
                        orig_conseq=os.path.join(folder, "{0}.{1}.conseq".format(record,alias))
                        for line in open(orig_conseq):
                                if first:
                                        if line[0]!=">": #Hack to add the description
                                                print(">{0}".format(record), file=conseq)
                                        first=False
                                print(line, end='', file=conseq)
                        commandline.append("-c={0}".format(conseq.name))
                        conseq.flush()
        if add_estseq:
                if not args.est_test_alias:
                        estseq=os.path.join(folder, "{0}.estseq.fa".format(record))
                else:
                        estseq=os.path.join(folder, "{0}.{1}.estseq.fa".format(record, args.est_test_alias))
                commandline.append("-e={0}".format(estseq))
        if ref_gtf!=None:
                commandline.append("-T={0}".format(ref_gtf))

        iscan=subprocess.Popen(commandline, shell=False, stdout=subprocess.PIPE)

        tempGtf=tempfile.NamedTemporaryFile(suffix=".gtf", mode='w')
        prefix=re.sub("\.gtf$", "", tempGtf.name)
        for line in iscan.stdout:
                print(line.decode('UTF-8').rstrip(), file=tempGtf)
        retcode=iscan.returncode
        # if retcode!=0:
        #       return out.name, retcode
        tempGtf.flush()
        subprocess.call([os.path.join(os.environ['EVAL_GTF'],'validate_gtf.pl'), '-f',
                         tempGtf.name, fasta], shell=False,
                        stdout=open(os.devnull,'w')) #, stderr=open(os.devnull,'w')) #I do not want any stupid output clogging the logs
        tempGtf.close() #Remove temporary file

        curr=None

        hits={}

        for line in open("{0}.fixed.gtf".format(prefix)):
                if line[0]=="#": continue
                fields=line.rstrip().split("\t")
                if fields==[] or fields=='' or len(fields)!=9: continue
                try:
                        if not re.search("codon", fields[2]): continue
                except IndexError:
                        raise IndexError(fields)
                tid=re.search(r".*transcript_id \"([^\"]*)\".*", fields[-1]).groups()[0]
                if tid not in hits:
                        hits[tid]={}
                        hits[tid]['start']=False
                        hits[tid]['stop']=False
                hits[tid][re.sub("_codon","",fields[2])]=True

        good_hits=[tid for tid in list(hits.keys()) if hits[tid]['start']==True and hits[tid]['stop']==True]
        
        for line in open("{0}.fixed.gtf".format(prefix)):
                if line[0]=="#":
                        print(line, end='', file=out)
                        continue
                fields=line.rstrip().split("\t")
                if len(fields)!=9: continue
                tid=re.search(r".*transcript_id \"([^\"]*)\".*", line).groups()[0]
                if tid in good_hits:
                        print(line, end='', file=out)
                # print(line, end='', file=out)

        
                #Check that the transcript has both stop and start codon

        out.close()
        os.remove("{0}.fixed.gtf".format(prefix))
        if add_conseq:
                conseq.close()

        return out.name, retcode

def writeXml( template, selected, out, split_directory=".", datasets=None, alias='', est_alias=None, masked=False, exclude_est=False, nscan=False):
        '''This function parses the XML instance file template (template variable) and inserts
        the various files in the appropriate positions.
        split_directory is the datastore directory of the split files.'''

    ##The following commented lines provide a brief example of what to look for in the XML template.
    # <annotation_files>
    # </annotation_files> 
    # <seq_files type="dna">
    # </seq_files>
    # <seq_files type="cons">
    # </seq_files>
    # <seq_files type="est">
    # </seq_files>

        out_folder=os.path.abspath(os.path.dirname(out))

        comb_fasta=open( os.path.join(out_folder, 'combined_genome.fa'),'w')
        comb_gtf=open( os.path.join(out_folder,  'combined_genome.gtf'),'w')
        print(">genome", file=comb_fasta)
        if nscan:
                comb_align=open(os.path.join(out_folder, 'combined_genome.align'),'w')
        else:
                comb_conseq=open(os.path.join(out_folder, 'combined_genome.conseq'),'w')

        if exclude_est: print("EST excluded", file=sys.stderr)
        else:
                comb_estseq=open(os.path.join(out_folder, 'combined_genome.estseq.fa'), 'w')
                print(">genome", file=comb_estseq)

        current_position=0
        if nscan:
                my_sequences={}
                first=True

        for name,folder in selected:
                if masked:
                        fasta=os.path.join(folder, "{0}.masked.fa".format(name))
                        if nscan: #{
                                if alias!='':
                                        align=os.path.join(folder,  "{0}.{1}.masked.align".format(name,alias) )
                                else:
                                        align=os.path.join(folder, "{0}.masked.align".format(name,alias) )

                                if not os.path.exists(align):
                                        align="{0}.{1}.align".format(name,alias)
                                        align=os.path.join(folder, align)
                                        assert os.path.exists(align)
                                        
                        else: #{
                                

                                if alias!='':
                                        conseq=os.path.join(folder,  "{0}.{1}.masked.conseq".format(name,alias) )
                                else:
                                        conseq=os.path.join(folder,  "{0}.{1}.conseq".format(name,alias) )
                                if not os.path.exists(conseq):
                                        conseq="{0}.{1}.conseq".format(name,alias)
                                        conseq=os.path.join( folder, conseq)
                                        assert os.path.exists(conseq), (name, alias, conseq)
                                #}
                else:
                        fasta=os.path.join(folder, "{0}.fa".format(name))
                        if nscan:
                                
                                align=os.path.join(  folder, "{0}.{1}.align".format(name,alias) )
                                assert os.path.exists(align)
                                for sequence in SeqIO.parse(open(align)):
                                        if not first:
                                                assert sequence.id in my_sequences, (name, sequence.id)
                                        else:
                                                my_sequences[sequence.id]=''
                                                first=False
                                        my_sequences[sequence.id]+=str(sequence.seq).rstrip() #Add the sequence
                        else:
                                conseq=os.path.join( folder, "{0}.{1}.conseq".format(name,alias) )
                                
                                assert os.path.exists(conseq), conseq

                fasta=SeqIO.read(open(fasta),'fasta')
                fasta_seq=re.sub(r"[a-z]","N", str(fasta.seq))
                if nscan:
                        for sequence in SeqIO.parse(open(align),'fasta'):
                                if not first:
                                        assert sequence.id in my_sequences, (name, sequence.id)
                                else:
                                        my_sequences[sequence.id]=''
                                        first=False
                                assert len(sequence.seq)==len(fasta_seq)
                                my_sequences[sequence.id]+=str(sequence.seq).rstrip() #Add the sequence
                                #}
                                
#               print(len(my_sequences), file=sys.stderr)
                print(fasta_seq, end='', file=comb_fasta)
                if not nscan:
                        conseq="".join([line.rstrip() for line in open(conseq)])

                        print(conseq, end='', file=comb_conseq)
                

#               gtf=open(os.path.join(folder, "{0}.gtf".format(name)))
                if not nscan:
                        gtf=open(os.path.join(folder, "{0}.gtf".format(name)))
                        
                else:
                        gtf=open(os.path.join(folder, "{0}.noutr.gtf".format(name)))

                for line in gtf:
                        line=line.rstrip().split("\t")
                        line[0]="genome"
                        line[3]=int(line[3])+current_position
                        line[4]=int(line[4])+current_position
                        print(*line, sep="\t", file=comb_gtf)

                if not exclude_est:
                        if not est_alias:
                                estseq=SeqIO.read(open(os.path.join( folder, '{0}.estseq.fa'.format(name))),'fasta')
                        else:
                                estseq=SeqIO.read(open(os.path.join( folder, "{0}.{1}.estseq.fa".format(name,est_alias) )),'fasta')
                        print(str(estseq.seq).rstrip(), file=comb_estseq, end='')
                current_position+=len(fasta)
        comb_fasta.close()
        comb_gtf.close()
        if not nscan:
                comb_conseq.close()
        else:
                for sequence in my_sequences:
                        print(">genome {0}".format(sequence), file=comb_align)
                        print(my_sequences[sequence], file=comb_align)
                comb_align.close()
        if not exclude_est: comb_estseq.close()

        out=open(out, 'w')
        for line in open(template):
                if '<!DOCTYPE iPE_instance SYSTEM "iPE_instance.dtd">' in line: #End of header
                        print(line, end='', file=out)
                        print("<!-- Datasets used for this instance: {0}-->".format(", ".join([str(i) for i in datasets]) ), file=out ) 
                elif '<annotation_files>' in line: #Print the GTF files
                        print(line, end='', file=out)
                        print(os.path.abspath(comb_gtf.name), file=out)

                elif '<seq_files type="dna">' in line: #Print the DNA files
                        print(line, end='', file=out)
                        print(os.path.abspath(comb_fasta.name), file=out)
                    
                elif '<seq_files type="cons">' in line:
                        print(line, end='', file=out)
                        print(os.path.abspath(comb_conseq.name), file=out)
                elif '<seq_files type="malign">' in line:
                        print(line, end='', file=out)
                        print(os.path.abspath(comb_align.name), file=out)
                                
                elif '<seq_files type="est">' in line:
                        print(line, file=out)
                        if exclude_est: continue
                        else: print(os.path.abspath(comb_estseq.name), file=out)
                else:
                        print(line, end='', file=out)

        return


def testTwinScan(outDir=".",test_dataset=[],args=None,pool=None, lock=None, debug=False):

        assert os.path.exists(outDir), outDir
                        ##output files
        outBare=os.path.join(os.getcwd(),outDir, "Bare")
        if not os.path.exists(outBare):
                os.makedirs(outBare)
        
        outConseq=os.path.join(os.getcwd(),outDir, "Conseq")
        if not os.path.exists(outConseq):
                os.makedirs(outConseq)
        
        outFull=os.path.join(os.getcwd(),outDir, "Full")
        if not os.path.exists(outFull):
                os.makedirs(outFull)

        zhmm=os.path.join(os.getcwd(), outDir, "parameters.zhmm")
        
        reference_list=open(os.path.join(os.getcwd(), outDir, "reference.list"),'w')
        jobs=dict()

        keys=['ref', 'bare', 'conseq', 'full']
        stripped_args=copy(args)
        del stripped_args.control_dataset
        del stripped_args.prepare_log
        print(stripped_args, file=sys.stderr)

        for name,folder in test_dataset:
                #                       tup_gtf=os.path.join(folder, "{0}.gtf".format(name))
                jobs[name]=dict().fromkeys(keys)

                if args.exclude_est:
                        tup_gtf=os.path.join(folder,"{0}.gtf".format(name))
                else:
                        tup_gtf=os.path.join(folder,"{0}.noutr.gtf".format(name))
                jobs[name]['ref']=tup_gtf

                if pool!=None:
#                        print(folder, name, zhmm, outBare, args, lock, file=sys.stderr)

                        jobs[name]['bare']=pool.apply_async(twinScan, 
                                                            args=(folder, name, zhmm, outBare, stripped_args, lock),
                                                            kwds={'add_conseq': False, 'add_estseq': False, 'ref_gtf': tup_gtf}) 
                        jobs[name]['conseq']=pool.apply_async(twinScan,
                                                              args=(folder, name, zhmm, outConseq, stripped_args, lock),
                                                              kwds={'add_conseq': True, 'add_estseq': False, 'ref_gtf': tup_gtf})
                        jobs[name]['full']=pool.apply_async(twinScan,
                                                            args=(folder, name, zhmm, outFull, stripped_args, lock),
                                                            kwds={'add_conseq': True, 'add_estseq': True, 'ref_gtf': tup_gtf})
                else:
                        jobs[name]['bare']=twinScan(folder, name, zhmm, outBare, args, lock, add_conseq=False, add_estseq=False, ref_gtf=tup_gtf)
                        jobs[name]['conseq']=twinScan(folder, name, zhmm, outConseq, args, lock, add_conseq=True, add_estseq=False, ref_gtf=tup_gtf)
                        jobs[name]['full']=twinScan(folder, name, zhmm, outFull, args, lock, add_conseq=True, add_estseq=True, ref_gtf=tup_gtf)

        #Print out
        bare_list=open(os.path.join(outDir, "{0}.bare.txt".format(args.test_prefix)),'w')
        conseq_list=open(os.path.join(outDir, "{0}.conseq.txt".format(args.test_prefix)),'w')
        full_list=open(os.path.join(outDir, "{0}.full.txt".format(args.test_prefix)),'w')
                       
        for name in jobs:
                if pool!=None:
                        for key in [key for key in keys if key!="ref"]:
#                                print(name, key)
                                jobs[name][key]=jobs[name][key].get()

                retcodes=[ jobs[name][key][1] for key in [key for key in keys if key!="ref"]]
                if [x for x in retcodes if x not in (0, None)]!=[]:
                        print("{0} not correctly evaluated; retcodes: {1}".format(name, retcodes), file=sys.stderr)
                        continue #Error in the evaluation!
                print(jobs[name]['ref'], file=reference_list)
                print(jobs[name]['bare'][0], file=bare_list)
                print(jobs[name]['conseq'][0], file=conseq_list)
                print(jobs[name]['full'][0], file=full_list)
                        
        bare_list.close()
        conseq_list.close()
        full_list.close()
        reference_list.close()

def trainTwinScan(outDir):
        '''This function calls ipestimate using the XML configuration file.'''
        os.chdir(outDir)
        if os.path.exists("parameters.zhmm"): return 0
        log=open("ipestimate.log",'w')
        retcode=subprocess.call(['ipestimate','instance.xml'], shell=False, stdout=log, stderr=log)
        return retcode

def argument_parser(argv):

        def check_negative(value):
                try: value=int(value)
                except:
                        raise argparse.ArgumentTypeError("{0} is not an integer!".format(value))
                if value<0:
                        raise argparse.ArgumentTypeError("{0} is not a positive integer ".format(value))
                return value

        def check_procs(value):
                try: value=int(value)
                except:
                        raise argparse.ArgumentTypeError("{0} is not an integer!".format(value))
                if value<0:
                        raise argparse.ArgumentTypeError("{0} is not a positive integer".format(value))
                if value>cpu_count():
                        raise argparse.ArgumentTypeError(
                                "This machine has at most {0} processors, you asked for {1}".format(cpu_count(), value))
                return value


        defaults={ 'twinscan': "/opt/N-SCAN/",
                   'test_prefix': 'predictions.dataset',
                   'nData': 4,
                   'proc': 4,
                   }

        parser=argparse.ArgumentParser(description="This script automates TwinScan training. It performs the splitting into k different datasets . The k datasets are then partitions in k instances, for each k-1 sets will be used for training and the k-th set will be used for testing. The test dataset will give name to the output folder.")
        testing=parser.add_mutually_exclusive_group()
        testing.add_argument('--control_dataset', type=argparse.FileType('r'), default=None, help="Optional prepare_data log with the selected models for testing.")
        testing.add_argument('-n', '--nData',  type=check_negative, default=10, help="Order of the K-cross validation. Default: 10" )
        parser.add_argument('--debug', action="store_true", default=False, help="Debug flag. If activated, multiprocessing will be disabled.")
        parser.add_argument('--wd', type=str, default="Instances", help="The folder where the instances will be saved. Default: '.' ")
        parser.add_argument('--split_log', type=str, required=True, help="The file into which the assignments are written into. Default: STDOUT. If a name is given, the file will be written into the working directory.")
        parser.add_argument("-p", "--processors", type=check_procs, dest="procs", default=defaults['proc'], help="Number of processors to use.")
        parser.add_argument('--masked', action="store_true", default=False, help="Flag. If set, *masked* sequences will be used for training.")
        parser.add_argument('--alias', type=str, required=True, help="The informant database alias.")
        parser.add_argument('--test_alias', type=str, help="[Optional] An alternative alias for CONseq files, for testing purposes. Defaults to alias.")
        parser.add_argument('--exclude_est', action="store_true", default=False, help="Flag. If set, ESTseq files will not be used.")
        parser.add_argument('--est_alias', type=str, default=None, help="[Optional] Alias for the ESTSeq files.")
        parser.add_argument('--est_test_alias', type=str, default=None,
                            help="[Optional] An alternative alias for ESTseq files, for testing purposes. Defaults to est_alias.")
        parser.add_argument('--twinscan', default=defaults['twinscan'], help="The TwinScan directory. Default: {0}".format(defaults['twinscan']))
        parser.add_argument('--nscan', default=False, action='store_true', help="Flag. If set, NScan (with .align files) will be performed.")
        parser.add_argument('--test_prefix', default=defaults['test_prefix'], help="The prefix name for the TwinScan prediction file.")
        parser.add_argument("--test_only", action="store_true", default=False, help="Flag. If called, the program will only try to perform the evaluation tests.")
        parser.add_argument("--dataset_master", default=None, help="Folder of prepare_dataset data. Default: same as the prepare_log file.")
        parser.add_argument('template', type=argparse.FileType('r'), help="The XML instance template for TwinScan.")
        parser.add_argument('--force', action="store_true", default=False)
        parser.add_argument('prepare_log', type=argparse.FileType('r'),
                            help="The dataset log of the prepare_datasets script. Used to infer the absolute paths of the files")
        args=parser.parse_args(argv)

        return args

#################

def main():
        args=argument_parser(sys.argv[1:])
        args.template.close()
        args.template=args.template.name
        if args.control_dataset!=None:
                args.nData=1
        prefix=os.path.join
        master=os.getcwd()
        if not os.path.exists(os.path.join(args.wd, "TwinScan")):
                os.makedirs(os.path.join(args.wd, "TwinScan"))
    
        with open(os.path.join(args.wd, 'twinscan_command.txt') , 'w') as command_out:
                print(*textwrap.wrap(
                                " ".join(sys.argv), break_on_hyphens=False, break_long_words=False ),
                       sep=" \\\n", file=command_out) #Nice printout

        if os.path.dirname(args.split_log)=="":
                args.split_log=os.path.join(args.wd, args.split_log)
        else: pass
    
        print("###Command line used:", file=sys.stderr)
        print(""," ".join(sys.argv), sep="\t", file=sys.stderr)

        ##This section creates the TWINSCAN environment variable, which is needed by twinscan.
        assert os.path.exists(args.twinscan) and os.path.isdir(args.twinscan) and os.path.exists("/".join([args.twinscan,"bin/iscan"]))

        os.environ['TWINSCAN']=args.twinscan
        #Add the TwinScan libraries to the PERL5 libraries
        if 'PERL5LIB' in os.environ:
                os.environ['PERL5LIB']+=":{twin}".format(twin=args.twinscan)
        else:
                os.environ['PERL5LIB']=":{twin}".format(twin=args.twinscan)
        args.twinscan=os.path.join(args.twinscan,"bin/iscan")
    
        ##If test_aliases are not defined they default to the training ones.
        if not args.test_alias:
                args.test_alias=args.alias[:]
                if args.masked: args.test_alias=args.test_alias #+".masked"

        if not args.est_test_alias and args.est_alias!=None: args.est_test_alias=args.est_alias[:]
        if args.est_test_alias=="None": args.est_test_alias=None
        
        if args.dataset_master==None:
                dataset_master = os.path.dirname(os.path.abspath(args.prepare_log.name)) #If just the file name is given, path is assumed to be the curent directory.
        else:
                args.dataset_master=os.path.abspath(args.dataset_master)
                assert os.path.exists(args.dataset_master) and os.path.isdir(args.dataset_master)
                dataset_master=args.dataset_master
        

        datasets=geneid.select(args.prepare_log, nData=args.nData, log=args.split_log, master=dataset_master) #Assign each gene model to a different dataset
        if args.control_dataset:
                control_master=os.path.dirname(os.path.abspath(args.control_dataset.name))
                control_dataset=geneid.select(args.control_dataset, nData=1,
                                              log=os.path.join(master, args.wd, "control_dataset.log"),
                                              master=control_master)        

        manager=Manager()
        lock=manager.Lock()
        
        if not args.debug:
                pool=Pool(processes=args.procs)
        else:
                pool=None
            
        jobs=[]
        for ds in range(args.nData):
                print(ds, file=sys.stderr)
                outDir=os.path.join(args.wd,"TwinScan","Test_{0}".format(ds))
                if not os.path.exists(outDir): os.makedirs(outDir)
                selected=set()
                selected_datasets=list(range(args.nData))
                if args.nData>1: selected_datasets.remove(ds)
                for key in selected_datasets:
                        selected=set.union(selected, datasets[key]) #Create the selected dataset
                selected=list(selected) #We need it to be ORDERED, not with random/rapid access.
                if args.force:
                        if os.path.exists(os.path.join(outDir, "parameters.zhmm")):
                                os.remove(os.path.join(outDir, "parameters.zhmm"))
                if not os.path.exists(os.path.join(outDir, "parameters.zhmm")):
                        outXml="/".join([outDir, "instance.xml"])
                        writeXml(args.template, selected, outXml, datasets=selected_datasets, 
                                 alias=args.alias, est_alias=args.est_alias, masked=args.masked, exclude_est=args.exclude_est, nscan=args.nscan)
                        if not args.debug:
                                job=pool.apply_async(trainTwinScan, args=(outDir,))
                                jobs.append([ds, job])
                        else:
                                print(trainTwinScan(outDir), file=sys.stderr)
                else: pass
        if not args.debug:
                jobs=[[job[0], job[1].get()] for job in jobs]

                print("Job", "Retcode", file=sys.stderr)

        for job in jobs: print(*job, file=sys.stderr)

        if args.nData>1 or args.control_dataset!=None:
                testJobs=[]
                print("##Begin performance assessment", file=sys.stderr)

                for ds in range(args.nData):
                        outDir=os.path.join(master, args.wd,"TwinScan" ,"Test_{0}".format(ds))
                        assert os.path.exists(outDir), outDir
                        if args.nData>1:
                                test_dataset=datasets[ds]
                        else:
                                test_dataset=control_dataset[0]

                        testTwinScan(outDir=outDir, test_dataset=test_dataset, args=args, pool=pool, lock=lock)
                        
                # for job in testJobs:
                #       if not args.debug:
                #               recordId, retcode = job.get()
                #       else:
                #               recordId, retcode = job
                # else: pass

                if pool!=None:
                        pool.close()
                        pool.join()

                

                #Do evaluate_gtf
                for ds in range(args.nData):

                        gtf_file=os.path.join(master, args.wd, "TwinScan", "Test_{0}".format(ds), "reference.list")
                        outBare=os.path.join(master, args.wd, "TwinScan", "Test_{0}".format(ds), "{0}.bare.txt".format(args.test_prefix))
                        outEval=open( os.path.join(master, args.wd, "TwinScan", "Test_{0}".format(ds), "eval.bare.txt"), 'w')
                        eval_gtf=subprocess.Popen([os.path.join(os.environ['EVAL_GTF'],'evaluate_gtf.pl'),'-Aq', gtf_file, outBare], stdout=outEval)
                        eval_gtf.communicate()
                        outEval.close()

                        outConseq=os.path.join(master, args.wd, "TwinScan", "Test_{0}".format(ds), "{0}.conseq.txt".format(args.test_prefix))
                        outConseqEval=open( os.path.join(master, args.wd, "TwinScan", "Test_{0}".format(ds), "eval.conseq.txt"), 'w')
                        eval_gtf=subprocess.Popen([os.path.join(os.environ['EVAL_GTF'],'evaluate_gtf.pl'),'-Aq', gtf_file, outConseq], stdout=outConseqEval)
                        eval_gtf.communicate()
                        outConseqEval.close()

                        if not args.exclude_est:
                                outFull=os.path.join(master, args.wd, "TwinScan", "Test_{0}".format(ds), "{0}.full.txt".format(args.test_prefix))
                                outFullEval=open( os.path.join(master, args.wd,"TwinScan", "Test_{0}".format(ds), "eval.full.txt"), 'w')
                                eval_gtf=subprocess.Popen([os.path.join(os.environ['EVAL_GTF'],'evaluate_gtf.pl'),'-Aq', gtf_file, outFull], stdout=outFullEval)
                                eval_gtf.communicate()
                                outFullEval.close()


                for word in ["bare", "conseq", "full"]:
                        outKcross=open(os.path.join(master, args.wd, "TwinScan", "kcross_mean.{0}.txt".format(word)),'w' )
                        evals=[os.path.join(os.getcwd(), args.wd, "TwinScan", "Test_{0}".format(i), "eval.{0}.txt".format(word)) for i in range(args.nData)]
                        for line in kcross(evals): print(line, file=outKcross)
                        outKcross.close()
   

if __name__=='__main__': main()
