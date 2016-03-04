#!/usr/bin/env python3


import sys
import subprocess
import os
import re
import shutil
from shutil import rmtree
import argparse
import textwrap


def select(split_log, master, out):
    '''This function retrieves the GB files (printed to the "out" file)
    and returns  the number of selected gene models.'''

    count = 0
    for line in split_log:
        name, folder, status = line.rstrip().split()[:3]
        if status!="FINISHED": continue
        count+=1
        folder=os.path.join(master, folder)
        for line in open(os.path.join(folder, name+".gb")):
            print(line, end='', file=out)
    return count

def main():

    def split_comma(string):
        return string.split(",")


    defaults={
              "augustus": "/opt/augustus.2.6.1/config/",
              "species": "eggplant",
              'proportion': 0.1,
              'out': 'models',
              'out_dir': 'Augustus',
              'aug_log': 'augustus_training',
              'kfold': 10,
              'opt_rounds': 4
        }


    parser=argparse.ArgumentParser()
    parser.add_argument("-a" ,"--augustus", type=str, default=defaults["augustus"], help="Augustus config path. Default: {0}".format(defaults['augustus']))
    parser.add_argument("--dataset_master", default=None, help="Folder of prepare_dataset data. Default: same as the prepare_log file.")
    parser.add_argument('--ests', type=argparse.FileType('r'), help="EST PSL file to use for training augustus." )
    parser.add_argument('--comparison', type=split_comma, default=[], help="Species to test the training against.")
    control=parser.add_mutually_exclusive_group()
    control.add_argument('--proportion', type=float, default=defaults['proportion'], help="Proportion of gene models to keep out for testing. Default: {0}".format(defaults['proportion']))
    control.add_argument("--control_dataset", type=argparse.FileType('r'), default=None, help="Optional prepare_data log with the selected models for testing.")
    aug_config=parser.add_argument_group("Options for augustus training.")
    aug_config.add_argument("--test_only", default=False, action="store_true", help="Flag. If set, the program will only evalaute the performances of <species> and eventually of other species provided by the \"comparison\" parameter.")
    aug_config.add_argument("-s", '--species', type=str, required=True, help="The species to train for. Required. It will be created starting from the tomato parameters.")
    aug_config.add_argument('--aug_log', type=str, default=defaults['aug_log'], help="Prefix for augustus training logs.")
    aug_config.add_argument('--kfold', type=int, default=defaults['kfold'], help="Number of validations for training Augustus. Equals the cpus used. Default: {0}".format(defaults['kfold']))
    aug_config.add_argument('--opt_rounds', type=int, default=defaults['opt_rounds'], help="Number of optimization rounds for Augustus. Default: {0}".format(defaults['opt_rounds']))
    parser.add_argument('--out', type=str, default=defaults['out'], help='The prefix of the output file. Default: {0}'.format(defaults['out']))
    parser.add_argument('--out_dir', type=str, default=defaults['out_dir'], help="Directory where the training will be performed. Default: {0}".format(defaults['out_dir']))
    parser.add_argument("prepare_log", type=argparse.FileType('r'), help="The log file from prepare_data.")
    args=parser.parse_args()

    if args.control_dataset!=None: args.proportion=0

    assert os.path.exists(args.augustus)
    assert os.path.exists(os.path.join(args.augustus,"species"))

    os.environ['AUGUSTUS_CONFIG_PATH']=args.augustus

    if args.proportion>=1 or (args.proportion==0 and not args.control_dataset):
        raise ValueError("The proportion value must be comprised between 0 and 1")

    master=os.getcwd()
    args.prepare_log=open(os.path.abspath(args.prepare_log.name))

    if args.dataset_master==None:
        dataset_master = os.path.dirname(os.path.abspath(args.prepare_log.name)) #If just the file name is given, path is assumed to be the curent directory.
    else:
        args.dataset_master=os.path.abspath(args.dataset_master)
        assert os.path.exists(args.dataset_master) and os.path.isdir(args.dataset_master)
        dataset_master=args.dataset_master


    if not os.path.exists(args.out_dir): os.makedirs(args.out_dir)
    os.chdir(args.out_dir)
    
    #Write the command line
    with open("augustus_command.txt",'w') as command_log:
        print(*textwrap.wrap(
                " ".join(sys.argv), break_long_words=False),
               sep="\\\n", file=command_log)

#    models_test_gtf=open(args.out+".test.gtf",'w')
    if not args.test_only:
        models=open(args.out+".gb",'w')
        num_models=select(args.prepare_log, dataset_master, models)
        for_test=int(round(num_models*args.proportion))
        models.close()

    if not args.control_dataset:

        split_proc=subprocess.Popen([os.path.join(
            args.augustus, "..",
            "scripts", "randomSplit.pl"), models.name, str(for_test)], shell=False)
        stdout, stderr=split_proc.communicate()
        if not split_proc.returncode==0:
            raise ValueError(stderr)
        models_test_list=list()
        for line in open(args.out+".gb.test"):
            if re.match("LOCUS", line):
                models_test_list.append(line.rstrip().split()[1].split("_")[0])

        for line in open(args.prepare_log.name):
            model, folder = line.rstrip().split()[:2]
            folder = os.path.join(os.path.dirname(args.prepare_log.name), folder)

    else:
        if not args.test_only:
            shutil.copyfile(args.out+".gb",args.out+".gb.train")
        model_test_gb=open(args.out+".gb.test", 'w')
        for line in args.control_dataset:
            model,folder,status=line.rstrip().split("\t")[:3]
            if status!="FINISHED": continue
            folder=os.path.abspath(os.path.join(
                    os.path.dirname(os.path.abspath(args.control_dataset.name)), folder))
            model_gb=os.path.join(folder, "{0}.gb".format(model))
            model_gtf=os.path.join(folder, "{0}.noutr.gtf".format(model))
            for line in open(model_gb): print(line.rstrip(), file=model_test_gb)

    if not args.test_only:

        #Remove previous version, if it exists
        if os.path.exists(os.path.join(args.augustus,"species", args.species)):
            rmtree(os.path.join(args.augustus,"species", args.species))

        if os.path.exists("autoAugTrain"): rmtree("autoAugTrain")
        
        commandline=[ os.path.join(args.augustus, "..", "scripts", "autoAugTrain.pl"),
                      "-v","-v","-v",
                      "--kfold={0}".format(args.kfold), "--cpus={0}".format(args.kfold),
                      "--CRF", "--optrounds={0}".format(args.opt_rounds),
                      "--trainingset={0}".format(models.name+".train"),
                      "--species={0}".format(args.species)]
        print(*textwrap.wrap(" ".join(commandline), break_long_words=False ), file=sys.stderr)
        if args.ests:
            commandline.append("--est={0}".format(args.ests))
    
        train_log=open(args.aug_log+".stdout.log", 'w')
        train_err=open(args.aug_log+".stderr.log", 'w')

        aug_train=subprocess.Popen(commandline,
                                   shell=False,
                                   stdout=train_log,
                                   stderr=train_err)

        aug_train.communicate()

      
    
    test_log=open(".".join([args.out, args.species, "txt"]), 'w')
    test_err=open(".".join([args.out, args.species, "stderr", "txt"]), 'w')

    
    test=subprocess.Popen([
            "augustus", "--species={0}".format(args.species),
            ".".join([args.out, "gb", "test"])], shell=False,
                          stdout=test_log,
                          stderr=test_err)

    test.communicate()
    test_log.close()
    test_gtf=open(".".join([args.out, args.species, "gtf"]),'w')
    models_test_gtf=open(args.out+".test.gtf",'w')
    for line in open(".".join([args.out, args.species, "txt"])):
        if line[0]=="#": continue
        if "a-posteriori" in line: break
        fields=line.rstrip().split("\t")
        if fields[1]=="database": print(line, end='', file=models_test_gtf)
        elif fields[2]=="AUGUSTUS": print(line, end='', file=test_gtf)
        
    test_gtf.close()
    models_test_gtf.close()

    eval_out=open(".".join([args.out, args.species, "eval"]),'w')
    eval_err=open(".".join([args.out, args.species, "eval.err"]),'w')

    eval_proc=subprocess.Popen([
            os.path.join(os.environ['EVAL_GTF'],"evaluate_gtf.pl"),
            '-Agq', args.out+".test.gtf",
            test_gtf.name], shell=False,
                               stdout=eval_out, stderr=eval_err)


    if args.comparison!=[]:
        for comparison in args.comparison:

            assert os.path.exists(
                os.path.join(args.augustus, "species", comparison)
                )
            comparison_log=open(".".join([args.out, comparison, "txt"]), 'w')
            comparison_err=open(".".join([args.out, comparison, "stderr", "txt"]), 'w')
            comparison_test=subprocess.Popen([
                    "augustus", "--species={0}".format(comparison),
                    ".".join([args.out,"gb","test"])], shell=False,
                                             stdout=comparison_log,
                                             stderr=comparison_err)
            comparison_test.communicate()
            comparison_log.close()
            comparison_gtf=open(".".join([args.out, comparison, "gtf"]), 'w')
            for line in open(".".join([args.out, comparison, "txt"])):
                if line[0]=="#": continue
                if "a-posteriori" in line: break
                fields=line.rstrip().split("\t")
                if fields[2]=="AUGUSTUS": print(line, end='', file=comparison_gtf)

            comparison_gtf.close()
            eval_out=open(".".join([args.out, comparison, "eval"]),'w')
            eval_err=open(".".join([args.out, comparison, "eval.err"]),'w')

            eval_proc=subprocess.Popen([
                    os.path.join(os.environ['EVAL_GTF'],"evaluate_gtf.pl"),
                    '-Agq', args.out+".test.gtf",
                    comparison_gtf.name], shell=False,
                                       stdout=eval_out, stderr=eval_err)


if __name__=='__main__': main()
