#!/usr/bin/env python

import sys, logging,re,os
import random
import io

def select(prepare_log, nData=0, log=sys.stdout, master="."):
    '''This function is used to create the datasets from the prepare_log file.'''

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
