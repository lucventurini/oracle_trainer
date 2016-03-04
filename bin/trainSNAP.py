#!/usr/bin/env python3


import sys
import re
import os
import subprocess
import argparse
import tempfile
from multiprocessing import Pool, Manager, cpu_count
from shutil import rmtree
import itertools
import textwrap
from Oracle_utils import geneid, zff2gtf
from Oracle_utils.kcross_mean import kcross


def trainSnap(records=None, args=None, workDir=None):
    '''This function automatizes the training for SNAP. The steps are as follows:
    1- recover and combine the .ann and fasta files from the record folder for each record.
    2- Use fathom on the combined ann and fasta files.
    3- Export the univoque and correct gene models
    4- Use forge to infer the parameters
    5- Use hmmassembler to construct the HMM file.'''

    #Command 0: create the GFF3 + sequence

    if not os.path.exists(workDir): os.makedirs(workDir)
    os.chdir(workDir)
    snap_hmm=os.path.join("{0}.hmm".format(args.model_name))
    if os.path.exists(snap_hmm):
        if not args.force: return snap_hmm

    paramFolder=os.path.join("snap_params")
    if os.path.exists(paramFolder): rmtree(paramFolder) #remove the params folder if it exists.
    os.makedirs(paramFolder)

    #Command 1: retrieve the .ann files present in the preparation folder

    annotFile=open("genome.ann", 'w')
    seqFile=open("genome.dna", 'w')

    for record in records:
        name, folder= record
        if args.masked:
            record_fasta=os.path.join(folder, "{0}.masked.fa".format(name))
        else:
            record_fasta=os.path.join(folder, "{0}.fa".format(name))
        for line in open(record_fasta): print(line, end='', file=seqFile)
        print(file=seqFile) #Print an empty line at the end of each sequence
        record_annot=os.path.join(folder, "{0}.ann".format(name))
        for line in open(record_annot): print(line, end='', file=annotFile)

    annotFile.close()
    seqFile.close()

    #Command 2: fathom

    fathom_log=open("fathom.log",'w')
    fathom=subprocess.Popen(['fathom', '-categorize', str(args.flank), annotFile.name, seqFile.name], 
                            shell=False, stdout=fathom_log, stderr=fathom_log)
    fathom.communicate()

    #Command 3: export
    export_log=open("export.log",'w')
    fathom_export=subprocess.Popen(['fathom', '-export', str(args.flank), '-plus', 'uni.ann', 'uni.dna'],
                                   shell=False, stdout=export_log, stderr=export_log)
    fathom_export.communicate()

    #Command 4: forge
    os.chdir(paramFolder)
    forge_log=open("forge.log",'w')
    forge=subprocess.Popen(['forge', '-lcmask', '../export.ann', '../export.dna'],
                           shell=False, stdout=forge_log, stderr=forge_log)
    forge.communicate()
    os.chdir("..")

    #Command 5: assemble
    snap_hmm=open(snap_hmm,'w')
    assemble_log = open("assemble.log", 'w')
    assemble=subprocess.Popen(['hmm-assembler', '-Z', args.family, args.model_name, 'snap_params'],
                              shell=False, stdout=snap_hmm, stderr=assemble_log)
    assemble.communicate()
    return snap_hmm.name

def testSnap(records=[], args=None, workDir=None, hmm=None ):
    
    '''This function has the purpose to check a SNAP parameter file
    and verify its precision against the test dataset.'''

    os.chdir(workDir)
    if not os.path.exists("snap_test"): os.makedirs("snap_test")
    ref_list=open("ref_list",'w')
    pred_list=open("pred_list",'w')

    for record in records:
        name, folder = record
        out=os.path.join("snap_test", "{0}.gtf".format(name))
        ref_gtf=os.path.join(folder, "{0}.noutr.gtf".format(name))
        fasta=os.path.join(folder, "{0}.masked.fa".format(name))
        if (not os.path.exists(out)) or os.stat(out).st_size==0:
            out=open(out,'w')
            tempGtf=tempfile.NamedTemporaryFile(suffix=".gtf")
            prefix=re.sub("\.gtf$", "", tempGtf.name)
            snap_comm=subprocess.Popen(['snap', '-lcmask', '-quiet', hmm, fasta ], 
                                       shell=False, stdout=subprocess.PIPE)
            
            
            for line in zff2gtf.zff2gtf(snap_comm.stdout):
                print(line, file=tempGtf)
            tempGtf.flush()
            subprocess.call(['/opt/eval-2.2.8/validate_gtf.pl', '-f',
                             tempGtf.name, fasta], shell=False,
                            stdout=open(os.devnull,'w'), stderr=open(os.devnull,'w')) #I do not want any stupid output clogging the logs
            tempGtf.close() #Remove temporary file

            for line in open("{0}.fixed.gtf".format(prefix)):
                print(line, end='', file=out)

            out.close()

        else: out=open(out)
        print(os.path.abspath(out.name), file=pred_list)
        print(os.path.abspath(ref_gtf), file=ref_list)
    ref_list.close()
    pred_list.close()

    evalFile=open("evaluation.txt",'w')
    eval_comm=subprocess.Popen(
        ['/opt/eval-2.2.8/evaluate_gtf.pl', '-qA', ref_list.name, pred_list.name],
        shell=False, stdout=evalFile)
    eval_comm.communicate()
    return
    #zffToGff3(zff) #I need to convert from *ZFF*. See http://www.bioperl.org/wiki/ZFF

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
        

    parser=argparse.ArgumentParser("This script will train automatically the data for SNAP.")
    parser.add_argument('-m', '--model_name', type=str, required=True, help="The name of the HMM model.")
    testing=parser.add_mutually_exclusive_group()
    testing.add_argument('--control_dataset', type=argparse.FileType('r'), default=None, help="Optional prepare_data log with the selected models for testing.")
    testing.add_argument('-n', '--nData',  type=check_negative, default=10, help="Order of the K-cross validation. Default: 10" )
    parser.add_argument('-p', '--processors', type=check_procs, default=0, help="Number of processors to use. Defaults to nData.")
    parser.add_argument('-f', '--family', type=str, choices=["worm","fly","plant"])
    parser.add_argument('--force', action='store_true', default=False,
                        help="Flag. If set, recalculation will take place.")
    parser.add_argument('--debug', action="store_true", default=False,
                        help="Flag. If active, it will disable multiprocessing to facilitate debugging.")
    parser.add_argument('--wd', required=True, help="Working directory.")
    parser.add_argument('--split_log', default="snap.log", help="Log to write the division to.")
    parser.add_argument('--masked', action="store_true", default=False,
                        help="Flag. If set, the maskedsequence will be used for training also.")
    parser.add_argument('--flank', type=int, default=2000, help="Flanking region to train SNAP. Default: 2000")
    parser.add_argument('--test_flank', type=int, default=2000, help="Flanking region to test SNAP. Default: 2000")
    parser.add_argument("--dataset_master", default=None, help="Folder of prepare_dataset data. Default: same as the prepare_log file.")
    parser.add_argument('prepare_log', type=argparse.FileType('r'), help="The preparation log from prepare_datasets")
    
    args=parser.parse_args(argv)

    return args

def main():

    args=argument_parser(sys.argv[1:])

    if args.control_dataset!=None:
        args.nData=1

    if args.processors==0:
        args.processors=min(args.nData,cpu_count())

    master=os.getcwd()
    if args.dataset_master==None:
        dataset_master = os.path.dirname(os.path.abspath(args.prepare_log.name)) #If just the file name is given, path is assumed to be the curent directory.
    else:
        args.dataset_master=os.path.abspath(args.dataset_master)
        assert os.path.exists(args.dataset_master) and os.path.isdir(args.dataset_master)
        dataset_master=args.dataset_master

    if not os.path.exists(args.wd):
        os.makedirs(args.wd)

    with open(os.path.join(args.wd, 'snap_command.txt') , 'w') as command_out:
        print(*textwrap.wrap(
                " ".join(sys.argv), break_long_words=False, break_on_hyphens=False),
               sep=" \\\n", file=command_out) #Nice printout

    if os.path.dirname(args.split_log)=="":
        args.split_log=os.path.join(master,args.wd,args.split_log)
        print(args.split_log, file=sys.stderr)

    datasets=geneid.select(args.prepare_log, nData=args.nData, log=args.split_log, master=dataset_master) #Assign each gene model to a different dataset
    if args.control_dataset:
        control_master=os.path.dirname(os.path.abspath(args.control_dataset.name))
        control_dataset=geneid.select(args.control_dataset, nData=1,
                                      log=os.path.join(master, args.wd, "control_dataset.log"),
                                      master=control_master)

    manager=Manager()
    lock=manager.RLock()
    pool=Pool(processes=args.processors)

    if not args.debug: jobs=[]
    hmms={}
    
    if args.masked:
        prefix=os.path.join(master, args.wd, "SNAP_masked")
    else:
        prefix=os.path.join(master, args.wd, "SNAP")

    for ds in range(args.nData):
        os.chdir(master)
        #Retrieve all the records which are NOT in the ds dataset
        if args.nData>1:
            records=list(set(
                    itertools.chain.from_iterable( 
                        datasets[x] for x in [y for y in range(args.nData) if y!=ds]
                        )))
        else:
            #Complete training
            records=datasets[ds]

        workDir=os.path.join(prefix, "Test_{0}".format(ds))

        if args.debug:
            snap_hmm=trainSnap(args=args, workDir=workDir, records=records)
            hmms[ds]=snap_hmm
        else:
            #def trainSnap(records, args, workDir=".")
            job=pool.apply_async(trainSnap, kwds={'workDir': workDir,
                                            'records': records,
                                            'args': args})
            jobs.append(tuple([ds,job]))

    if not args.debug:
        for job in jobs:
            ds, snap_hmm = job
            hmms[ds]=snap_hmm.get()
    
    if args.nData>1:
        #Testing
        if not args.debug: test_jobs=[]
        for ds in range(args.nData):
            os.chdir(master) #Necessary
            workDir=os.path.join(prefix, "Test_{0}".format(ds))

            test_set=datasets[ds]
            if args.debug:
                testSnap(records=test_set, args=args, workDir=workDir, hmm=hmms[ds] )
            else:
                job=pool.apply_async( testSnap, kwds={'records': test_set,
                                                      'args': args,
                                                      'workDir': workDir,
                                                      'hmm': hmms[ds]
                                                      })
                test_jobs.append(job)

        for job in test_jobs:
            job.get()
        outKcross=open(os.path.join(prefix, "kcross_mean.txt"),'w')
        evals=[os.path.join(prefix, "Test_{0}".format(ds), "evaluation.txt") for ds in range(args.nData)]
        for line in kcross(evals): print(line, file=outKcross)
        outKcross.close()
                           

    elif args.control_dataset:
        if not args.debug: test_jobs=[]
        print("Controlling", file=sys.stderr)
        os.chdir(master)
        workDir=os.path.join(prefix, "Test_0")
        
        test_set=control_dataset[0]
        print(len(test_set), file=sys.stderr)
        if args.debug:
            testSnap(records=test_set, args=args, workDir=workDir, hmm=hmms[0])
        else:
            print(hmms[0], len(test_set), workDir, args.__dict__, file=sys.stderr)
            job=pool.apply_async( testSnap, kwds={'records': test_set,
                                                  'args': args,
                                                  'workDir': workDir,
                                                  'hmm': hmms[0]
                                                  })

    if not args.debug:
        for job in test_jobs:
            job.get()

        pool.close()
        pool.join()


if __name__=='__main__': main()
