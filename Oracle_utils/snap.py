#!/usr/bin/env python

from Oracle_utils.interface import Interface
import sys,subprocess,os,logging
import operator
from Bio import SeqIO
import textwrap

class SNAP(Interface):

    def __init__(self, interface=None):
        if interface is None:
            raise ValueError("No interface provided!")

        self.inherit(interface) #Inherit logger, datafolder, etc.
        self.define_step()

    def define_step(self):
        self.step="SNAP"

    def __call__(self):
        self.prepareZff()
        super(self.__class__, self).__call__()


    def prepareZff(self):
        '''This method creates a ZOE file as required by SNAP.'''

        self.logger.info("Starting the ZFF creation.")
        
        gff3=os.path.join(self.dataFolder, "{0}.gff3".format(self.new_name))

        gff3_lines = [line.rstrip() for line in open(gff3) ]
        gff3_lines = list(filter( lambda line: line[0]!="#" and line.split()[2]=="CDS",
                             gff3_lines))

        if len(gff3_lines)==0:
            self.logger.error("No GFF3 lines! Exiting.")
            self.failed=True
            return

        zff = os.path.join(self.dataFolder, "{0}.ann".format(self.new_name))
        if os.path.exists(zff) and os.path.getsize(zff)>0:
            self.logger.warn("ZFF file already present. Exiting.")
            return

        self.logger.debug("Creating and sorting the ZFF lines.")
        lines = []
        strand=None

        for line in gff3_lines:
            line=line.rstrip().split("\t")
            if not strand: strand = line[6]
            line[3]=int(line[3])
            line[4]=int(line[4])
            lines.append(["Exon", line[3], line[4], self.new_name])


        try:
            if strand!="-":
                self.logger.debug("Plus strand. Normal sorting.")
            #Plus strand. Order first by start than end
                lines=sorted(
                    sorted(lines, key=operator.itemgetter(2)),
                    key=operator.itemgetter(1))

            else:
                self.logger.debug("Minus strand. Reverse sorting, and reverse fields.")
            #Minus strand. Order first by end than start
                new_lines = []
                for line in lines:
                    new_lines.append([line[0], line[2], line[1], self.new_name])
                lines=new_lines[:]
                lines = sorted(lines, key=operator.itemgetter(1), reverse=True)
                if len(lines)==1:
                    lines[0][0]="Esngl"
                else:
                    lines[0][0]="Einit"
                    lines[-1][0]="Eterm"
        except Exception as error:
            self.logger.error("Problems in sorting!")
            self.logger.exception(error)
            self.failed=True
            return


        zff_file=open(zff, 'w')
        print(">{0}".format(self.new_name), file=zff_file)
        for line in lines:
            print(*line, sep="\t", file=zff_file)
        zff_file.close()
        self.logger.info("Finished creating the ZFF")
        return

        
