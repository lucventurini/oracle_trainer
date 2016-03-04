#!/usr/bin/env python3

import subprocess
import argparse
import hashlib
from Bio import SeqIO
from multiprocessing import Pool, Process, Manager
import logging
import logging.handlers  # Need *explicit* import
import textwrap
import shutil
from Oracle_utils import augustus, twinscan, snap, geneid, interface, initializer
from Oracle_utils.utils import *  # Import the file creation functions
from collections import OrderedDict as odict


#######################################

def manageSeq(record, args, new_name, genomeLengths, main_queue):
    '''This is the master function. It takes the CD-HIT sequence and produces all necessary files.'''

    # Initialize the process
    name=record.id
    folderName=hashlib.md5(name.encode("UTF-8")).hexdigest().upper()
    folderName=os.path.join(*[args.master]+splitString(folderName, n=2)[:2]+[name])
    dataFolder=os.path.join(os.getcwd(),folderName)

    #Attacch to the main logger
    main_logger = logging.getLogger('main_{name}'.format(name=name))
    main_logger.propagate = False #Otherwise messages will be sent to STDERR also
    main_logger.setLevel(args.debug)
    queue_handler = logging.handlers.QueueHandler(main_queue)
    main_logger.addHandler(queue_handler) #The default settings will take care of what will be written or not to the main log.
    
    main_logger.info("Started with {0}".format(new_name))

    if not os.path.exists(dataFolder):
        main_logger.debug("Created folder {0} for {1}".format(dataFolder, new_name))
        os.makedirs(dataFolder)
    elif os.path.exists(dataFolder) and args.force:
        main_logger.debug("Removed folder {0} for {1}".format(dataFolder, new_name))
        shutil.rmtree(dataFolder)
        os.makedirs(dataFolder)

    #Create the local logger
    local_logger=logging.getLogger("log_{name}".format(name=name))
    local_handler=logging.FileHandler(os.path.join(dataFolder, "log.txt"))
    formatter=logging.Formatter("%(asctime)s - %(module)s - %(levelname)s ---- %(message)s")
    local_handler.setFormatter(formatter)
    if args.debug == logging.DEBUG:
        local_handler.setLevel(args.debug)
    else:
        local_handler.setLevel(logging.INFO)

    local_logger.addHandler(local_handler)
    local_logger.setLevel(args.debug)
    local_logger.propagate=False

    status="FINISHED"

    #Create the Interface

    current_interface = interface.Interface(logger=local_logger,
                                            name = name,
                                            dataFolder = dataFolder,
                                            new_name = new_name,
                                            record=record,
                                            namespace=args)


    current_initializer = initializer.Initializer(
        interface = current_interface)
    current_initializer()

    if current_initializer.failed:
        status="FAILED"

    module_dict = odict(
        [
        ("snap", snap.SNAP),
        ("augustus", augustus.Augustus),
        ("geneid", geneid.GeneID),
        ("twinscan", twinscan.TwinScan)
        ]
        )

    for mod in module_dict:
        #As soon as something goes wrong, break.
        if status=="FAILED":
            break
        current_mod = module_dict[mod](interface = current_initializer)
        current_mod()

        if current_mod.failed==True:
            status="FAILED"

    if status!="FAILED":
        local_logger.info("Finished preparation.")
        main_logger.info("Finished {0}".format(new_name))
    else:
        main_logger.error("{0} failed at the {1} level".format(
                new_name, current_initializer.step))

    del local_logger

    return new_name, folderName, status

    

############################################


def getWuBlastParams(blast_param):

    '''This function parses the WUBLAST parameter file.
    The said file must contain one flag per line;
    at least one of them must contain the alias information.
    E.g.:
         alias   tomato
    Flags must be preceded by a dash, e.g.:
         -lcmask'''


    params={'options': list()}
    for line in blast_param:
        if line[0]=="-":
            params['options'].append(line.rstrip().split()[0])
            continue
        line=line.rstrip().split("\t")
        if line[0]=="alias": params['alias']=line[1]
        elif line[0]=="database":
            params['database']=line[1]
            if not os.path.exists(params['database']+".xnt"): subprocess.call(['/opt/wublast/xdformat', '-n', '-I', params['database']], shell=False)
        else:
            params['options'].append("=".join(line[:2]))

    assert 'alias' in params, "An alias must be provided for WuBlast in the parameter file!"

#    params['options']=" ".join(params['options'])

    return params

############################################

def getMegaBlastParams(blast_param):


    params={'options': list()}
    for line in blast_param:
        # if line[0]=="-":
        #     params['options'].append(line.rstrip().split()[0])
        #     continue
        #Avoid line formatting breaking everything
        line=[x for x in line.rstrip().split("\t") if x!=''] 
        if line[0]=="alias": params['alias']=line[1]
        elif line[0]=='database':
            params['options']+=['-db',line[1]]
            params['database']=line[1]
            if not os.path.exists(line[1]+".nsq"):
                subprocess.call(['makeblastdb', '-in', line[1],
                                 '-parse_seqids', '-dbtype', 'nucl'],
                                shell=False, stderr=subprocess.PIPE )
        elif len(line)==1:
            params['options'].append(line[0])
        else:
            params['options']+=["{0}".format(line[0]), line[1]]
                                 
    return params

############################################

def argument_parser(argv):

    '''This function parses the command line and restitutes the proper arguments in an argparse oject.'''
    defaults={ 'nData': 4,
               'master': 'datasets',
               'p': 10,
               'flank': 2000,
               'log': sys.stdout,
               'start_from': 1,
               }

    parser=argparse.ArgumentParser('This program prepares all the data necessary for the training of the various gene predictors.')
    parser.add_argument('-f','--genomic_fasta', type=argparse.FileType('r'), required=True, help="The FASTA file of the scaffolds.")
    parser.add_argument('-m', '--masked_fasta', type=argparse.FileType('r'), required=True, help="The *masked* FASTA file of the scaffolds. Used for testing.")
    parser.add_argument('-p','--processors', type=int, default=defaults['p'], help="Number of processors to use. Default: {0}".format(defaults['p']))
    parser.add_argument('-o', '--result_file', type=argparse.FileType('a'), default=defaults['log'], help="The log file on which all results will be written. Default: stdout")
    parser.add_argument('-l', '--log_file',  default=sys.stderr,
                        help="Log file. Default: STDERR")

    parser.add_argument("--force", action="store_true", default=False,
                        help="Flag. If set, all calculations will be repeated from scratch.")

    parser.add_argument("--only_canonical", action="store_true", default=False,
                        help="Flag; if set, only transcripts with canonical junctions will be retained.")

    parser.add_argument('--twinscan', type=str, default="/opt/N-SCAN/", help="The TWINSCAN directory.")
    parser.add_argument('-g','--genomic_gtf', required=True, dest='gengtf', type=argparse.FileType('r'), help='The GTF file of the annotations.')
    parser.add_argument('--wublast', default=None, type=argparse.FileType('r'), help="WU-Blast parameter file.")
    parser.add_argument('--megablast', default=None, type=argparse.FileType('r'), help="NCBI-Blast parameter file.")
    parser.add_argument('--last', default=None, type=argparse.FileType('r'), help="LAST genome file. **ATTENTION**: the FASTA file should be a single (not multi) FASTA with the name equating the name of the sequence. E.g. \"arabidopsis\" with a single sequence, \">arabidopsis\". ")
    parser.add_argument('--conseq', default=None, type=argparse.FileType('r'), help="Pre-generated Conseq file.")

    parser.add_argument('--conseq_alias', default=None, type=str, help="Alias for the Conseq files.")

    # blat=parser.add_mutually_exclusive_group()
    parser.add_argument('--blat_target', default=None, type=argparse.FileType('r'), help="BLAT database for alignment.")
    parser.add_argument('--estseq', default=None, type=argparse.FileType('r'), help="Pre-generated EstSeq file.")
    parser.add_argument('--estseq_alias', default=None, type=str, help="Alias to use for estseq files when importing from an external file.")

    parser.add_argument('--flank', type=int, default=defaults['flank'], help="The amount of genomic sequence to flank the genomic loci. Default: {0}".format(defaults['flank']))
    parser.add_argument('--padding', type=int, default=60, help="Amount of sequence to chunk up for GeneID signals.")
    parser.add_argument('-d', '--datastore', dest='master', type=str, default=defaults['master'], help="The master folder for all the datasets. Default: {0}".format(defaults['master']))
    parser.add_argument('--debug', choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="WARNING", help="Debug verbosity level.")
    parser.add_argument('--geneid', action='store_true', default=False, help="Prepare the files necessary for GeneID.")
    parser.add_argument('--start_from', default=defaults['start_from'], type=int, help="Number from which the counting for models will begin from. Useful for adding gene models to an existing dataset." )
    parser.add_argument('cdhit', type=argparse.FileType('r'), help="The CD-Hit clustered file.")

    args=parser.parse_args(argv)

    args.debug=getattr(logging, args.debug)

    return args


def queue_listener_process(main_queue=None):
    '''Simple listener to the queue. It logs the messages to the main logger, and exits when the right signal
    END QUEUE NOW
    is provided.'''
    logger=logging.getLogger('main')

    while True:
        record=main_queue.get()
        if record=="END QUEUE NOW": break
        logger.handle(record)


def main():

    args=argument_parser(sys.argv[1:]) #Parse the arguments

    #Create the logger
    logger=logging.getLogger("main")

    formatter=logging.Formatter("%(asctime)s - %(module)s - %(levelname)s ---- %(message)s")

    if args.log_file!=sys.stderr:
        handler=logging.FileHandler(args.log_file, mode='w')
        # stream_handler=logging.StreamHandler()
        # stream_handler.setFormatter(formatter)
        # stream_handler.setLevel(logging.ERROR)
        # logger.addHandler(stream_handler)
    else:
        handler=logging.StreamHandler()

    handler.setFormatter(formatter)
    handler.setLevel(args.debug)
    logger.addHandler(handler)
    logger.setLevel(args.debug)
    logger.propagate=False


    with open('prepare_dataset.command' , 'w') as command_out:
        print(*textwrap.wrap(
                " ".join(sys.argv), break_on_hyphens=False, break_long_words=False ),
               sep=" \\\n", file=command_out) #Nice printout


    logger.info("Loading the previous dataset.")
    previous_dataset={}
    previous_dataset['models']=set()
    if os.path.exists(args.result_file.name):
        for line in open(args.result_file.name):
            model_name, folder = line.rstrip().split("\t")[:2]
            name=os.path.basename(folder)
            previous_dataset[name]=[model_name]
            previous_dataset['models'].add(model_name)
        args.result_file=open(args.result_file.name, 'w')
            

    manager=Manager()
    lock=manager.RLock()
    main_queue=manager.Queue(-1)

    #First part: let's start with the cdhit FASTA file...
    logger.debug("Checking the TWINSCAN configuration.")
    if (not os.path.exists(args.twinscan)) or \
            (not os.path.isdir(args.twinscan)) or \
            (not os.path.exists( os.path.join(args.twinscan, "bin", "iscan"))):
        logger.critical("TwinScan not found!")
        sys.exit(1)
    os.environ['TWINSCAN']=args.twinscan

    if not os.path.exists(args.master):
        os.makedirs(args.master)


    logger.debug("Creating the reference files.")
    args.genomic_fasta.close()
    args.genomic_fasta=args.genomic_fasta.name
    args.masked_fasta.close()
    args.masked_fasta=args.masked_fasta.name
    if args.last!=None:
        args.last.close()
        args.last=args.last.name

    genomeLengths=dict()
    for fasta_file in (args.genomic_fasta, args.masked_fasta):
        if ( not os.path.exists(fasta_file+".db") ) or \
                os.stat(fasta_file+".db").st_ctime<os.stat(fasta_file).st_ctime:
            logger.debug("Creating the DB for the reference genome.")
            db=SeqIO.index_db(fasta_file+".db", filenames=[fasta_file], format='fasta')
            if genomeLengths==dict():
                genomeLengths=dict((chrom, len(db[chrom])) for chrom in db)
            
    setattr(args, "genomic_fasta", "{0}.db".format(args.genomic_fasta))
    setattr(args, "masked_fasta", "{0}.db".format(args.masked_fasta))


    logger.debug("Loading BLAST parameters.")
    if args.wublast:
        args.wublast=getWuBlastParams(args.wublast)
    if args.megablast:
        args.megablast=getMegaBlastParams(args.megablast)

    if args.conseq:
        if not args.conseq_alias:
            logger.critical("An alias is required if you are loading a pre-calculated Conseq file!")
            sys.exit(1)
        args.conseq.close()
        SeqIO.index_db("{0}.db".format(args.conseq.name), 
                       filenames=[args.conseq.name], format="fasta")
        args.conseq="{0}.db".format(args.conseq.name)

    if args.wublast and args.conseq:
        if args.wublast['alias']==args.conseq_alias:
            logger.critical("You cannot provide the same alias for the internal wublast and the external CONSeq!")
            sys.exit(1)
        if "{0}.masked".format(args.wublast['alias'])==args.conseq_alias:
            logger.critical("You cannot provide the same alias for the internal wublast and the external CONSeq!")

    if args.blat_target:
        args.blat_target.close()
        args.blat_target=args.blat_target.name
    if args.estseq:
        if not args.estseq_alias:
            logger.critical("An alias for external import of ESTSeq files is needed!")
            sys.exit(1)
        args.estseq.close()
        logger.debug("Indexing the external ESTSeq file.")
        SeqIO.index_db("{0}.db".format(args.estseq.name), 
                       filenames=[args.estseq.name], format="fasta")
        args.estseq="{0}.db".format(args.estseq.name)

        
    args.gengtf.close()
    args.gengtf=args.gengtf.name

    logger.info("Starting the preparation of data.")

    jobs=[]
    count=args.start_from-1

    pool=Pool(processes=args.processors) #Create the pool of processors

    #Queue listener, to log the relative messages to the defined log file.
    queue_listener = Process(target=queue_listener_process,
                             kwargs = {'main_queue': main_queue})
    queue_listener.start()
        
    stripped_args=strip_namespace(args) #From oracle_utils

    logger.warn("Starting the batch analysis.")
    for record in SeqIO.parse(args.cdhit, 'fasta'):
        if record.id in previous_dataset:
            new_name=previous_dataset[record.id][0]
        else:
            count+=1
            new_name="MODEL{0}".format(count)
            while new_name in previous_dataset['models']:
                count+=1
                new_name="MODEL{0}".format(count)
            previous_dataset[record.id]=[new_name]
            previous_dataset['models'].add(new_name)

        job=pool.apply_async(manageSeq, args=(record,stripped_args,new_name, genomeLengths, main_queue))
        jobs.append(job)

    for job in jobs:
        print(*job.get(), sep="\t", file=args.result_file)

    main_queue.put_nowait("END QUEUE NOW")

    pool.close() #Close the pool
    pool.join()

    queue_listener.join() #Close the listener

    logger.warn("Run finished.")
        

if __name__=='__main__': main()
