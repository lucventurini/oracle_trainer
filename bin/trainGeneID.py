#!/usr/bin/env python3


import sys
import os
import subprocess
import random
import textwrap
import shutil
from multiprocessing import Pool, Manager, cpu_count
import argparse
from math import ceil
from Oracle_utils import geneid_utils
from Oracle_utils.kcross_mean import kcross
from copy import copy

def createTraining(records, args, folder, lock):
    '''This function will call the prepare_template() and optimize() functions from the geneid module. It will return the putatively optimized version of the parameter file.'''

#    print(records, file=sys.stderr)
    #prepare_template(records, args, lock, out)

    unoptim=os.path.join(folder, "unoptimized.param")
    if args.clean:
        if os.path.exists(unoptim): os.remove(unoptim)
        #Remove previous calculations
        for filename in [filename for filename in os.listdir(folder) if "_freq" in filename or "_entropy" in filename]:
            os.remove(os.path.join(folder,filename))


    if not os.path.exists(unoptim) or os.stat(unoptim).st_size==0 or args.force:
        unoptim=open(unoptim, 'w')
        for line in geneid_utils.prepare_template(records, folder, args):
            print(line, file=unoptim) #Fill in the values
        unoptim.close()
    else:
        unoptim=open(unoptim)
        unoptim.close()

def optimizeTraining(records, args, folder, lock, pool=None):

    unoptim=open(os.path.join(folder, "unoptimized.param")) #Bug in the code - I considered it as a file, with a "name" attribute
    optimized=os.path.join(folder, "optimized.param")
    if args.clean:
        if os.path.exists(os.path.join(folder, "GID_Optimization")):
            shutil.rmtree(os.path.join(folder, "GID_Optimization"))

    if args.force or args.clean:
        if os.path.exists(optimized): os.remove(optimized)

    print(os.path.join(folder, "GID_Optimization"), file=sys.stderr)
    if not os.path.exists(os.path.join(folder, "GID_Optimization")):
        print("Non-existent folder")
        os.makedirs(os.path.join(folder, "GID_Optimization"))

    if (not os.path.exists(optimized)) or (os.stat(optimized).st_size==0):
        optimized=open(optimized,'w')

        #Override runtime command 
        population_file=os.path.join(folder, "GID_Optimization", "population_records.log")

        if os.path.exists(population_file):
            optimization_records=set()
            with open(population_file) as pop:
                for line in pop:
                    name, folder = line.rstrip().split()
                    optimization_records.add( tuple([name, folder]) )
            lock.acquire()
            print("Previous population recovered, size: {0}".format(len(optimization_records)), file=sys.stderr)
            lock.release()
        else:
            if args.optimization_proportion<1 and args.optimization_proportion>0:
                pop_size=int(ceil(len(records)*args.optimization_proportion))
            elif args.optimization_size>0:
                assert args.optimization_size < len(records)
                print("Input size:", args.optimization_size, file=sys.stderr)
                pop_size=args.optimization_size
            else:
                pop_size=len(records)

            assert pop_size>0
            if len(records)==pop_size:
                optimization_records=records
            else:
                optimization_records=random.sample(records, pop_size)
            lock.acquire()
            print(pop_size, file=sys.stderr)
            lock.release()

            with open(population_file,'w') as out:
                for record in optimization_records:
                    print(*record, file=out)
        #####
        #Start optimization
        for line in geneid_utils.optimize(optimization_records, args, unoptim, lock=lock, pool=pool ):
            print(line, file=optimized) # fill in the values
    lock.acquire()
    print("Finished optimization for folder {0}".format(folder), file=sys.stderr)
    lock.release()

    return

def testTraining(dsFolder, records, args, param=None):
    '''This function tests the final parameter file.'''

    if not param:
        param=os.path.join(dsFolder, "optimized.param") #optimized.param
        assert os.path.exists(param), os.listdir(dsFolder)
        gtf_folder=os.path.join(dsFolder, "test_optimized")
        evalFile=os.path.join(dsFolder, ".".join(["test", "optimized", "eval"]))
        predList=os.path.join(dsFolder, ".".join(["test","optimized","list"]))
    else:
        gtf_folder=os.path.join(dsFolder, "test_{0}".format(os.path.basename(param)))
        evalFile=os.path.join(dsFolder, ".".join(["test", os.path.basename(param), "eval"]))
        predList=os.path.join(dsFolder, ".".join(["test", os.path.basename(param), "list"]))

    evalFile=open(evalFile,'w')
    refList=open(os.path.join(dsFolder, "test.reference.txt"  ),'w')
    predList=open(predList,'w')

    if os.path.exists(gtf_folder):
        from shutil import rmtree
        rmtree(gtf_folder)

    os.makedirs(gtf_folder)

    for record in records:
        name, folder=record 
        gtf=os.path.join(gtf_folder, "{0}.noutr.gtf".format(name))
        if not os.path.exists(gtf):
            fasta=os.path.join(folder,"{0}.masked.fa".format(name))
            assert os.path.exists(fasta)
            found=geneid_utils.doGeneID(gtf, param, fasta) #Has GID found any models for the sequence?
        else: pass #gtf=open(gtf)

        print(os.path.join(folder, name+".gtf"), file=refList)
        print(gtf, file=predList)

    refList.close()
    predList.close()              
    evaluation=subprocess.Popen([os.path.join(os.environ['EVAL_GTF'], 'evaluate_gtf.pl'), '-Aq',
                                refList.name, predList.name], stdout=evalFile, stderr=open(os.devnull,'w'))

    evaluation.communicate()
    evalFile.close()
    return geneid_utils.get_eval(evalFile.name)


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


    def maximum_one(string):
        try: val=float(string)
        except:
            raise ValueError("Pval must be a number between 0 and 1!")
        if val<=0 or val>1: raise ValueError("This must be a number between 0 and 1!")
        return val

    def positive(string):
        try:
            val=int(string)
            if val<3: raise ValueError
            return val
        except ValueError:
            raise ValueError("The CDS MM order must be a positive number higher than 3, you provided {0}".format(string))

    def percentage(string):
        try:
            val=float(string)
            if val<0:
                raise ValueError("This value must be greater than 0%!")
            if val>1:
                val/=100
            
            if val>1:
                raise ValueError("This value must be lower than 100%")
            return val
        except:
            raise


    defaults={ 'iEwf': -4.5,
               'fEwf': 1.5,
               'dEwf': 0.5,
               'iOwf': 0.2,
               'dOwf': 0.05,
               'fOwf': 0.45,
               'min_zscore': 0.2,
               'cds_maximum': 5,
               'signal_maximum': 1,
               'start_maximum': 2,
               'window_size': 5,
               'size_exon_diff': 0,
               'signal_cutoff': -5,
               'split_log': 'split.log',
               'optimization_size': 0,
               'optimization_proportion': 1,
        }

    parser=argparse.ArgumentParser("This is a full GeneID training pipeline.")
    parser.add_argument('-p','--processors', dest='proc', default=1, type=int, help="The number of processors to use.")
    parser.add_argument('-w', '--wd', type=str, help="The destination folder")
    parser.add_argument('--comparison_only', action='store_true', default=False, help="Flag. If set, the program will only perform the comparisons with available parameter files.")
    testing=parser.add_mutually_exclusive_group()
    testing.add_argument('--control_dataset', type=argparse.FileType('r'), default=None, help="Optional prepare_data log with the selected models for testing.")
    testing.add_argument('-n', '--nData',  type=check_negative, default=10, help="Order of the K-cross validation. Default: 10" )


    markov_estimation=parser.add_argument_group("Parameters for the creation of the unoptimized parameter file.")
    markov_estimation.add_argument('--signal_maximum', choices=[0,1,2], type=int, dest="signal_maximum",
                                   default=defaults['signal_maximum'],
                                   help="Maximum order of the signal MM. Default: {0}".format(defaults['signal_maximum']) )
    markov_estimation.add_argument('--start_maximum', choices=[0,1,2], type=int, dest="start_maximum",
                                   default=defaults['start_maximum'],
                                   help="Maximum order of the signal MM. Default: {0}".format(defaults['start_maximum']) )
    markov_estimation.add_argument('--cds_maximum', type=positive, dest="cds_maximum",
                                   default=defaults['cds_maximum'],
                                   help="Maximum order of the CDS/Intron MM. Must be higher than 2. Default: {0}".format(defaults['cds_maximum']) )


    markov_estimation.add_argument("--window_size", default=defaults['window_size'], type=check_negative,
                                   help="Size of the overlapping windows for the Z-score calculation. Default: {0}".format(defaults['window_size']))
    markov_estimation.add_argument("--min_zscore", default=defaults['min_zscore'], type=float,
                                   help="Minimum averaged z-score to consider a position as informative. Default: {0}".format(defaults['min_zscore']))
    markov_estimation.add_argument("--size_exon_diff", type=float, default=defaults['size_exon_diff'],
                                   help="Difference between size factor and exon factor. In tomato and potato, this is equal to 0.2. Default: {0}".format(
            defaults['size_exon_diff']))

    markov_estimation.add_argument("--start_cutoff", type=float, default=defaults['signal_cutoff'], 
                                   help="Minimum cutoff for considering a start signal. Default: {0}".format(defaults['signal_cutoff']))
    markov_estimation.add_argument("--acceptor_cutoff", type=float, default=defaults['signal_cutoff'], 
                                   help="Minimum cutoff for considering an acceptor signal. Default: {0}".format(defaults['signal_cutoff']))
    markov_estimation.add_argument("--donor_cutoff", type=float, default=defaults['signal_cutoff'], 
                                   help="Minimum cutoff for considering a donor signal. Default: {0}".format(defaults['signal_cutoff']))


    optimus=parser.add_argument_group('Optimization parameters.')
    optimus.add_argument('--no_optimization', default=False, action='store_true', help="Flag. If set, no optimization will be performed. Useful for testing the raw parameter file.")
    optimus.add_argument('--iEwf', default=defaults['iEwf'], type=float, help="Initial exon weight factor for testing. Default: {0}".format(defaults['iEwf']))
    optimus.add_argument('--dEwf', default=defaults['dEwf'], type=float, help="Exon weight factor step. Default: {0}".format(defaults['dEwf']))
    optimus.add_argument('--fEwf', default=defaults['fEwf'], type=float, help="Final exon weight factor for testing. Default: {0}".format(defaults['fEwf']))
    optimus.add_argument('--iOwf', default=defaults['iOwf'], type=float, help="Initial exon factor for testing. Default: {0}".format(defaults['iOwf']))
    optimus.add_argument('--dOwf', default=defaults['dOwf'], type=float, help="Exon factor step. Default: {0}".format(defaults['dOwf']))
    optimus.add_argument('--fOwf', default=defaults['fOwf'], type=float, help="Final exon factor for testing. Default: {0}".format(defaults['fOwf']))
    optimus.add_argument("--optimization_size", default=defaults['optimization_size'], type=check_negative,
                         help="Fixed size of the dataset for optimization. Useful when the training set is huge (>2-300 sequences). Default: None (the whole input dataset will be used for optimization).")
    optimus.add_argument("--optimization_proportion", default=defaults['optimization_proportion'], type=percentage,
                         help="Relative size of the dataset for optimization. Useful when the training set is huge (>2-300 sequences). Default: 100%% (the whole input dataset will be used for optimization).")

    parser.add_argument('--debug', action="store_true", default=False, 
                        help="If active, additional logging is performed and operations are NOT performed through multiprocessing.")
    parser.add_argument('--force', default=False, action="store_true", help="Flag. If set, everything will be recalculated.")
    parser.add_argument("--clean", default=False, action="store_true", help="Flag. If set, all files (except for FASTA files) will be cleaned in order to start from scratch.")
    parser.add_argument("--comparison", type=str, action='append', help="The reference parameter file(s) to compare against. Can be specified multiple times.")
    parser.add_argument('--masked', action='store_true', default=False, help="Use masked sequences to create the training.")
    parser.add_argument('--dataset_master', default=None, help="Folder of prepare_dataset data. Default: same as the prepare_log file.")
    parser.add_argument("--split_log", type=str, default=defaults['split_log'], help="The splitting log from trainTwinScan (if available).")
    parser.add_argument('template', type=argparse.FileType('r'), help="The template file")
    parser.add_argument('prepare_log', type=argparse.FileType('r'), help="The preparation log from prepare_dataset")
#    parser.add_argument('geneid', type=str, nargs='?', default=sys.stdout , help="The output file")

    if "-h" in argv or "--help" in argv:
        parser.print_help()
        parser.exit()
        return

    return parser.parse_args(argv)


def main():

    args=argument_parser(sys.argv[1:])
    if args.control_dataset!=None:
        args.nData=1

    args.template.close()
    args.template=args.template.name
    if not args.comparison: args.comparison=[]

    # if args.geneid!=sys.stdout:
    #     if not os.path.exists(os.path.dirname(os.path.abspath(args.geneid))):
    #         os.makedirs(os.path.dirname(os.path.abspath(args.geneid)))
    #     args.geneid=open(os.path.abspath(args.geneid), 'w')

    print("###Command line used:", file=sys.stderr)
    print(""," ".join(sys.argv), sep="\t", file=sys.stderr)
    if not os.path.exists(args.wd): os.makedirs(args.wd)
    master=os.getcwd()
    with open(os.path.join(args.wd, 'geneid_command.txt') , 'w') as command_out:
        print(*textwrap.wrap(
                " ".join(sys.argv), break_on_hyphens=False, break_long_words=False),
               sep=" \\\n", file=command_out) #Nice printout
        # print(""," ".join(sys.argv), sep="\t", file=command_out) #Print out the command used

    if os.path.dirname(args.split_log)=="":
        args.split_log=os.path.join(args.wd, args.split_log) #Add the working directory to the split_log


    if args.dataset_master==None:
        dataset_master = os.path.dirname(os.path.abspath(args.prepare_log.name)) #If just the file name is give, path is assumed to be the curent directory.
    else:
        args.dataset_master=os.path.abspath(args.dataset_master)
        assert os.path.exists(args.dataset_master) and os.path.isdir(args.dataset_master)
        dataset_master=args.dataset_master

    datasets=geneid_utils.select(args.prepare_log, nData=args.nData, log=args.split_log, master=dataset_master) #Assign each gene model to a different dataset
    if args.control_dataset:
        control_master=os.path.dirname(os.path.abspath(args.control_dataset.name))
        control_dataset=geneid_utils.select(args.control_dataset, nData=1,
                                      log=os.path.join(master, args.wd, "control_dataset.log"),
                                      master=control_master)        

#    if ("--processors" not in sys.argv or "-p" not in sys.argv) and args.nData>1: args.proc=args.nData

    manager=Manager()
    pool=Pool(processes=args.proc)
    lock=manager.RLock()

    backup_args=copy(args)
    print(args, file=sys.stderr)
    del args.prepare_log
    del args.control_dataset

    if not args.comparison_only:
        jobs=[]
        dataset_inputs=dict()

        for ds in range(args.nData):
            print(ds, file=sys.stderr)
            outDir=os.path.join(args.wd, "GeneID", "Test_{0}".format(ds))
            if not os.path.exists(outDir): os.makedirs(outDir)
            selected=set()
            selected_datasets=list(range(args.nData))
            if args.nData>1: selected_datasets.remove(ds)
            for key in selected_datasets:
                selected=set.union(selected, datasets[key]) #Create the selected dataset
            selected=list(selected) #We need it to be ORDERED, not with random/rapid access.
            dataset_inputs[ds]={'selected': selected, 'outDir': outDir}
            if not args.debug:
#                stripped_args=copy(args)
                jobs.append(pool.apply_async(createTraining, args=(selected, args, outDir, lock),
                                             ))

        if not args.debug:
            for job in jobs:
                job.get()

        if not args.no_optimization:#Start the optimization
            for ds in range(args.nData):
                optimizeTraining( dataset_inputs[ds]['selected'], args,
                                  dataset_inputs[ds]['outDir'], lock,
                                  pool=pool)

            
#            createTraining(selected, args, outDir, lock, pool=pool)

        

        # if not args.debug:
        #     for job in jobs: job.get()

    test_jobs=[]
    if args.nData>1:
        for ds in range(args.nData):
            folder=os.path.join(args.wd,"GeneID", "Test_{0}".format(ds))
            if args.debug:
                print(ds,   *testTraining(folder, datasets[ds], args, param=os.path.join(folder,"unoptimized.param")), sep="\t")
                if not args.no_optimization:
                    print(ds,   *testTraining(folder, datasets[ds], args), sep="\t")
                if args.comparison!=[]:
                    for param in args.comparison:
                        print(ds, param, *testTraining(folder, datasets[ds], args, param=param), sep="\t")
            else:
                test_jobs.append(tuple([ds, "unoptimized", pool.apply_async(testTraining, args=(folder, datasets[ds], args ),
                                                                            kwds={'param': os.path.join( folder, "unoptimized.param") }) ]) )

                if not args.no_optimization:
                    test_jobs.append(tuple([ds, "optimized", pool.apply_async(testTraining, args=(folder, datasets[ds], args ) )]) )
                
                if args.comparison!=[]:
                    for param in args.comparison:
                        test_jobs.append(tuple([ds, param, pool.apply_async(testTraining, args=(folder, datasets[ds], args), kwds={'param': param})]))


    elif backup_args.control_dataset!=None:
        folder=os.path.join(args.wd, "GeneID", "Test_0")
        test_set=control_dataset[0]
        ds=0
        if args.debug:
            print(ds, *testTraining(folder, test_set, args, param=os.path.join(folder, "unoptimized.param")), sep="\t")
            if not args.no_optimization:
                    print(ds,   *testTraining(folder, test_set, args), sep="\t")
            if args.comparison!=[]:
                for param in args.comparison:
                    print(ds, param, *testTraining(folder, test_set, args, param=param), sep="\t")
        else:
            test_jobs.append(tuple([ds, "unoptimized", pool.apply_async(testTraining, args=(folder, test_set, args ),
                                                                        kwds={'param': os.path.join( folder, "unoptimized.param") }) ]) )

            if not args.no_optimization:
                test_jobs.append(tuple([ds, "optimized", pool.apply_async(testTraining, args=(folder, test_set, args ) )]) )
                
            if args.comparison!=[]:
                for param in args.comparison:
                    test_jobs.append(tuple([ds, param, pool.apply_async(testTraining, args=(folder, test_set, args), kwds={'param': param})]))
                

    for job in test_jobs:
         print(job[0], job[1], *job[2].get())

    if backup_args.control_dataset or args.nData>1:
        outKcross=open(os.path.join(master, args.wd, "GeneID", "kcross_mean.unoptimized.txt"),'w')
        evals=[os.path.join(os.getcwd(), args.wd, "GeneID", "Test_{0}".format(i), "test.unoptimized.param.eval") for i in range(args.nData)]
        for line in kcross(evals): print(line, file=outKcross)
        outKcross.close()

    if not args.no_optimization:
        if backup_args.control_dataset or args.nData>1:
            outKcross=open(os.path.join(master, args.wd, "GeneID", "kcross_mean.optimized.txt"),'w')
            evals=[os.path.join(os.getcwd(), args.wd, "GeneID", "Test_{0}".format(i), "test.optimized.eval") for i in range(args.nData)]
            for line in kcross(evals): print(line, file=outKcross)
            outKcross.close()

    if args.comparison!=[]:
        for param in args.comparison:
            name=os.path.basename(param)
            outKcross=open(os.path.join(master, args.wd, "GeneID", "kcross_mean.{0}.txt".format(name)),'w')
            evals=[os.path.join(os.getcwd(), args.wd, "GeneID", "Test_{0}".format(i), "test.{0}.eval".format(name)) for i in range(args.nData)]
            for line in kcross(evals): print(line, file=outKcross)
            outKcross.close()

    pool.close()
    pool.join()


if __name__=='__main__': main()
