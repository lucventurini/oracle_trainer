#!/usr/bin/env python

from Oracle_utils.interface import Interface
import subprocess
import os
import textwrap


class Augustus(Interface):

    def __init__(self, interface=None):

        if interface is None:
            raise ValueError("No interface provided!")

        self.inherit(interface) #Inherit logger, datafolder, etc.
        self.define_step()


    def define_step(self):
        self.step="AUGUSTUS"
        

    def prepareGb(self, masked=False):
        '''Main function. It will create the GeneBank file. Two versions:
        masked, with the masked sequence, and vanilla.'''

        gff3=os.path.join(
                self.dataFolder, "{0}.gff3".format(self.new_name))

        if masked:
            self.logger.info("Creating masked GeneBank.")
            gb = os.path.join( self.dataFolder, "{0}.masked.gb".format(self.new_name))
            fasta=os.path.join(
                self.dataFolder, "{0}.masked.fa".format(self.new_name))


        else:
            self.logger.info("Creating un-masked GeneBank.")
            gb = os.path.join( self.dataFolder, "{0}.gb".format(self.new_name))
            fasta=os.path.join(
                self.dataFolder, "{0}.fa".format(self.new_name))

        command_line = [ os.path.join(
                os.environ['AUGUSTUS_CONFIG_PATH'], "..", "scripts", "gff2gbSmallDNA.pl"),
                         gff3,
                         fasta,
                         str(self.flank),
                         gb]
        
        gb_conversion = subprocess.Popen( command_line,
                                          shell=False,
                                          stdout=subprocess.DEVNULL,
                                          stderr=subprocess.PIPE)
        stdout, stderr = gb_conversion.communicate()
        if gb_conversion.returncode not in (0,None):
            stderr=stderr.decode("UTF-8")
            self.failed=True
            self.logger.critical("GeneBank conversion failed. STDerr:")
            self.logger.critical("\n".join(textwrap.wrap(
                        stderr, break_on_hyphens=False, break_long_words=False)
                                           ))
            return


    def __call__(self):
        '''Execution function. It will create a GB for both masked and unmasked sequences.'''

        for flag in (False, True):
            self.prepareGb(masked=flag)
            if self.failed==True: break

        super(self.__class__, self).__call__()
