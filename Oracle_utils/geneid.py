#!/usr/bin/env python

from Oracle_utils.interface import Interface
import sys,re,os
from myRecords.gtfRecord import gtfRecord
from myRecords import GTF
from Bio import SeqIO

class GeneID(Interface):

    def __init__(self, interface=None):
        
        if interface is None:
            raise ValueError("No interface provided!")

        self.inherit(interface) #Inherit logger, datafolder, etc.
        self.define_step()

    def define_step(self):
        self.step="GeneID"


    def __call__(self):
        for flag in [False,True]:
            self.prepareGeneId(masked=flag)
            if self.failed==True:
                break
        #Call the Interface call function to finish.
        super(self.__class__, self).__call__()


    def prepareGeneId(self, masked=False):

        if masked:

            self.logger.info("Starting the creation of GeneID files for the masked sequence.")
            addition=".masked"
            fasta=os.path.join( self.dataFolder, "{0}.masked.fa".format(self.new_name))


        else:
            addition="" #This is a string snippet used for formatting file names.
            fasta = os.path.join( self.dataFolder, "{0}.fa".format(self.new_name))
            self.logger.info("Starting the creation of GeneID files for the unmasked sequence.")

        fasta=SeqIO.read(open(fasta),'fasta')
        gtf = os.path.join( self.dataFolder, "{0}.gtf".format(self.new_name))
        self.logger.info("Analysing the GTF record.")
        current_record = gtfRecord(logger=self.logger)
        for gtfLine in GTF.GTF(gtf): current_record.add(gtfLine)

        #I have to import the values here.
        current_record.calculateSignals(fasta, length=self.padding,
                                        flank=10**6, strict=self.only_canonical)

        if current_record.failed:
            #Something went wrong. Exit.
            self.failed=True
            return

        self.logger.info("Finished analysing the GTF record, printing the various features.")


        for feature in "introns false_starts false_donors false_acceptors".split():
            self.logger.debug("Printing feature: {0}".format(feature))
            outFile=open(os.path.join(
                    self.dataFolder,
                    "{0}.{1}{2}.fa".format(self.new_name,feature, addition)),'w')
            for fastaFeature in current_record.__dict__[feature]:
                print(fastaFeature.format('fasta'), file=outFile)
            
            outFile.close()
    
        for feature in "donors acceptors".split():
            self.logger.debug("Printing feature: {0}".format(feature))
            outFile=open(os.path.join(
                    self.dataFolder,
                    "{0}.{1}{2}.fa".format(self.new_name,feature,addition)),'w')
            
            for fastaFeature in current_record.__dict__[feature]:
                print(
                    current_record.__dict__[feature][fastaFeature].format('fasta'),
                    file=outFile, end='')
            
            outFile.close()


        cdsFile=open(os.path.join(
                self.dataFolder, "{0}.cds{1}.fa".format(self.new_name, addition)),'w')
        self.logger.debug("Printing CDS")

        print(current_record.cds.format('fasta'), file=cdsFile, end='')
        cdsFile.close()
        self.logger.debug("Printing Stop sequence")
        stopFile=open(os.path.join(
                self.dataFolder, "{0}.stop{1}.fa".format(self.new_name, addition)),'w')
        print(current_record.stop_fasta.format('fasta'), file=stopFile, end='')

        self.logger.debug("Printing Start sequence")
        startFile=open(os.path.join(self.dataFolder, 
                                    "{0}.start{1}.fa".format(self.new_name, addition)),'w')
        print(current_record.start_fasta.format('fasta'), file=startFile, end='')
        self.logger.info("Finished the GeneID routine.")
        startFile.close()

        return
