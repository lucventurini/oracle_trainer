#!/usr/bin/env python

import sys,re,os
from copy import copy
import io


'''Generic utilities for the suite.'''


def splitString(string, n=2):
    '''Recursive function to split a sequence into chunks of a given length.
    It accepts the lenght of the fragments as keyword argument (n)'''

    if n<=0 or not isinstance(n, int):
        raise ValueError("The minimum length must be an integer of size greater than 0.")
    if len(string)<=n: return [string]
    else: return [string[:n]]+splitString(string[n:], n=n)


def select(prepare_log, nData=0, log=sys.stdout, master="."):
    '''This function is used to create or import the datasets for Kcross validation.'''

    if nData<=0: raise ValueError("nData must be greater than 0, you provided {0}!".format(nData))

    data=set()
    if not isinstance(prepare_log, io.IOBase):
        assert os.path.exists(prepare_log)
        prepare_log=open(prepare_log)

    for line in prepare_log:
        name, folder, status = line.rstrip().split()[:3]
        if status!="FINISHED": continue
        folder=os.path.join(master, folder)
        data.add(tuple([name, folder]))

    repopulate=False
    datasets=dict().fromkeys(list(range(nData)))
    for key in datasets: datasets[key]=set()

    if log!=sys.stdout and os.path.exists(log):
        for line in open(log):
            ds, name, folder = line.rstrip().split()
            ds=int(ds)
            if ds not in datasets:
                repopulate=True
                break
            datasets[ds].add(tuple([name, folder]))

        return datasets

    for val in data:
        dataset=random.sample(list(range(nData)),1)[0]
        datasets[dataset].add(val)

    if log!=sys.stdout:
        log=open(log, 'w')
    for ds in range(nData):
        for val in datasets[ds]:
            name, folder = val
            print(ds, name, folder, sep="\t", file=log)

    if log!=sys.stdout: log.close()
    return datasets

def strip_namespace(args):
    '''Quick utility to strip the ARGS NameSpace of all buffer-like handles. Necessary for multiprocessing.'''
    new_args=copy(args)
    for key in new_args.__dict__:
        if type(getattr(new_args,key)) is io.TextIOWrapper:
            setattr(new_args,key, getattr(new_args,key).name)
    return new_args


