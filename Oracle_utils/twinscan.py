#!/usr/bin/env python

from Oracle_utils.interface import Interface
import sys,re,os
from myRecords.gtfRecord import gtfRecord
from myRecords import GTF
from Bio import SeqIO
import subprocess
from Bio.Blast.NCBIXML import parse as xparser
from . import conseq
import gzip

class TwinScan(Interface):

    def __init__(self, interface=None):

        if interface is None:
            raise ValueError("No interface provided!")

        self.inherit(interface) #Inherit logger, datafolder, etc.
        self.define_step()

    def define_step(self):
        self.step="TwinScan"

    def __call__(self):

        grocery = {
            'conseq': [self.conseq, self.extract_conseq],
            'estseq': [self.estseq, self.extract_estseq],
            'megablast': [self.megablast, self.doMegaBlast],
            'wublast': [self.wublast, self.doWuBlast],
            }

        for action in grocery:
            if grocery[action][0] is not None:
                grocery[action][1]()
            if self.failed==True:
                break

        #Call the super calling method to close everythin.
        super(self.__class__, self).__call__()
        return
            

#########################
    def extract_conseq(self):
        '''Method to extract the CONSeq slice from a CONSeq genome file.'''


        self.logger.debug("Extracting the slice from CONSeq.")
        conseqOut = os.path.join(
            self.dataFolder, "{0}.{1}.conseq".format(self.new_name,
                                                     self.conseq_alias))
        if os.path.exists(conseqOut) and os.path.getsize(conseqOut)>0:
            self.logger.warn("CONSeq file for alias {0} already present. Exiting".format(self.conseq_alias))
            return

        try:
            index = SeqIO.index_db(self.conseq)
        except Exception as error:
            self.logger.exception(error)
            self.logger.error("Conseq is: {0}".format(self.conseq))
            self.failed=True
            return


        if self.chrom not in index:
            self.logger.critical("Chromosome {0} not in Conseq file. Exiting.".format(self.chrom))
            self.failed=True
            return

        #Establish boundaries
        start = max(self.start-self.flank-1, 0)
        end = min(len(index[self.chrom])-1, 
                  self.end-1+self.flank)
    
        seq=index[self.chrom][start:end]
        with open(conseqOut,'w') as out:
            print(seq.seq, file=out)

        self.logger.debug("Finished extracting the CONSeq sequence.")
        return

###########################

    def extract_estseq(self):

        '''Function to extract and print the ESTSeq slice.'''

        self.logger.debug("Extracting the slice from ESTSeq.")
        estseqOut = os.path.join(
            self.dataFolder, "{0}.{1}.fa".format(self.new_name,
                                                 self.estseq_alias))

        if os.path.exists(estseqOut) and os.path.getsize(estseqOut)>0:
            self.logger.warn("ESTSeq file for alias {0} already present. Exiting".format(
                    self.estseq_alias))
            return
    
        index=SeqIO.index_db(self.estseq)
        if self.chrom not in index:
            self.logger.critical("Chromosome {0} not in ESTSeq. Exiting.".format(self.chrom))
            self.failed=True
            return
        start = max(self.start-self.flank-1, 0)
        end = min(len(index[self.chrom])-1, 
                  self.end-1+self.flank)
    
        seq=index[self.chrom][start:end]
        with open(estseqOut,'w') as out:
            seq.id=self.new_name
            seq.description=""
            print(seq.format('fasta'), file=out)

        self.logger.debug("Finished extracting the CONSeq sequence.")
        return

###########################

    def doMegaBlast(self, masked=True):
        '''This function performs the WU-BLASTN and the conseq calling.
        Ideally it should be used to map the **masked** sequence against a **masked** genome.
        The method uses the custom conseq module to convert the MegaBLAST into a CONSeq file.'''

        if masked:
            self.logger.info("Starting masked MegaBLAST")
            fasta=os.path.join(self.dataFolder, "{0}.masked.fa".format(self.new_name))
            addition="masked."
        else:
            self.logger.info("Starting unmasked MegaBLAST")
            fasta=os.path.join(self.dataFolder, "{0}.fa".format(self.new_name))
            addition=''

        conseqOut=os.path.join(
            self.dataFolder,
            "{0}.{1}.{2}conseq".format(self.new_name, self.megablast['alias'], addition))

        if os.path.exists(conseqOut):
            self.logger.warn("WuBlast already performed. Exiting.")
            return True # avoid redoing the wuBlast

        blastout=os.path.join(
            self.dataFolder,
            "{0}.{1}.{2}megablast.xml.gz".format(
                self.new_name, self.megablast['alias'], addition))

        if not os.path.exists(blastout):
            blastout_file=gzip.open(blastout,'wt') #Open the blastout
            self.logger.info("Performing the MegaBLAST")
            command_line = ["blastn", "-query", fasta,
                            '-outfmt', '5']
            command_line += self.megablast['options']
            megablast_process = subprocess.Popen(command_line,
                                                 shell=False,
                                                 stdout=subprocess.PIPE,
                                                 stderr=subprocess.PIPE)
            stdout, stderr = megablast_process.communicate()
            if megablast_process.returncode not in (0,None):
                stderr=stderr.decode("UTF-8")
                os.remove(blastout_file)
                self.logger.critical("Error in the MegaBLAST.")
                self.logger.critical("STDERR: {0}".format(stderr))
                self.failed=True
                return
            stdout=stdout.decode("UTF-8")
            print(stdout, file=blastout_file)
            blastout_file.close()


        self.logger.debug("Loading the MegaBLAST XML")
        blast_record = next(xparser(gzip.open(blastout, 'rt')))
        conseqOut = open(conseqOut, 'w')
        self.logger.debug("Printing the CONSeq")
        print(conseq.conseq(
                blast_record, unaligned=True,
                gaps=False, transitions=False,
                with_name=False), file=conseqOut)
        return
        

###########################

    def doWuBlast(self, masked=True):

        ###Not implemented yet
        '''This function performs the WU-BLASTN and the conseq calling.
        It maps the **masked** sequence against a **masked** genome.'''

        if masked:
            self.logger.info("Starting masked WUBLAST")
            fasta=os.path.join(self.dataFolder, "{0}.masked.fa".format(self.new_name))
            addition="masked."
        else:
            self.logger.info("Starting unmasked WUBLAST".format(addition=addition))
            fasta=os.path.join(self.dataFolder, "{0}.fa".format(self.new_name))
            addition=''

        conseqOut=os.path.join(
            self.dataFolder,
            "{0}.{1}.{2}conseq".format(self.new_name, self.wublast['alias'], addition))

        if os.path.exists(conseqOut):
            self.logger.warn("WuBlast already performed. Exiting.")
            return True # avoid redoing the wuBlast

        blastout=os.path.join(
            self.dataFolder,
            "{0}.{1}.{2}blastn".format(
                self.new_name, self.wublast['alias'], addition))

        if os.path.exists(blastout+".gz"):
            self.logger.debug("Decompressing previously obtained WUBLAST.")
            subprocess.call(['gzip', '-d', blastout+".gz"], shell=False)
        elif not os.path.exists(blastout):
            self.logger.debug("Starting WUBLAST against {database}".format(database=self.wublast['database']))
            blastout=open(blastout,'w')
            command_line = [ os.path.join(os.environ['WUBLAST'], 'blastn'),
                                          self.wublast['database'], fasta]

            command_line += self.wublast['options']


            wublastProc=subprocess.Popen(command_line, shell=False, stdout=blastout, stderr=subprocess.PIPE)
            stdout, stderr = wublastProc.communicate()
            if wublastProc.returncode!=0:
                logger.error("WUBLAST failed")
                logger.error("STDERR:")
                stderr=stderr.decode("UTF-8")
                stderr="\n".join(textwrap.wrap(stderr, break_long_words=False, break_on_hyphens=False))
                logger.error(stderr)
                self.failed = True
                return

            blastout.close()
            blastout=blastout.name

        conseqOut=open(conseqOut,'w')

        self.logger.debug("Converting the WUBLAST into a CONSeq file.")
        consProcess=subprocess.Popen([os.path.join(os.environ['TWINSCAN'], 'bin', 'conseq.pl'),
                                      "-u", fasta, blastout],
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False) #Create the conseq

    ##Having a definition line breaks the *parameter estimation*.
    ##However, said definition line is necessary for the **prediction**
    ##and will have to be added at runtime.
        stdout, stderr = consProcess.communicate()
        if consProcess.returncode not in (0,None):
            stderr=stderr.decode("UTF-8")
            logger.error("Conversion to CONSeq failed")
            logger.error("STDERR: {0}".format(stderr))
            self.failed=True
            return
        stdout=stdout.decode("UTF-8").split("\n")
        for line in stdout:
            line=line.rstrip()
            if line=="": continue
            if line[0]==">": continue #print(">{0}".format(wublast['alias']))
            print(line, file=conseqOut, end='')

        self.logger.debug("Compressing the WUBLAST output.")
        comprBlast=subprocess.call(['gzip', blastout], shell=False)

#############################################

def doBlat(self, masked=True):

    '''This function performs BLAT and the creation of the ESTSeq file.
    As the arguments suggest, it maps against the masked file.'''

    fasta=os.path.join(
        self.dataFolder, "{0}.masked.fa".format(self.new_name))

    self.logger.info("Starting BLAT")

    estOut=os.path.join(
        self.dataFolder, "{0}.estseq.fa".format(self.new_name))
    blatout=os.path.join(
        self.dataFolder, "{0}.psl".format(self.new_name))

    if os.path.exists(estOut) and os.path.exists(blatout+".gz"):
        self.logger.warn("ESTSeq already present. Exiting.")
        return # avoid redoing the wuBlast

    estseqOut=open(estOut,'w')
    if os.path.exists(blatout+".gz"):
        self.logger.debug("Decompressing BLAT")
        subprocess.call(['gzip','-d', blatout+".gz"], shell=False)
    else:
        self.logger.debug("Performing BLAT")
        command_line = ['blat',
                        fasta,
                        self.blat_target,
                        blatout]
                        

        blat_command=subprocess.Popen(command_line,
                                      shell=False,
                                      stdout=subprocess.DEVNULL,
                                      stderr=subprocess.PIPE)
        stdout,stderr=blat.communicate()
        if blat.returncode not in (0, None):
            stderr=stderr.decode("UTF-8")
            stderr="\n".join(textwrap.wrap(
                    stderr, break_on_hyphens=False, break_long_words=False))

            logger.error("BLAT failed!")
            logger.error("STDERR: {0}".format(stderr))
            self.failed = True
            return

    
    self.logger.debug("Converting to ESTSeq.")
    command_line = [os.path.join(os.environ['TWINSCAN'],'bin', 'estseq.pl'),
                    masked_fasta, blatout]

    estseqProcess=subprocess.Popen(command_line,
                                   shell=False,
                                   stdout=estseqOut, stderr=subprocess.PIPE)

    stdout, stderr = estseqProcess.communicate()
    if estseqProcess.returncode not in (0,None):
        stderr=stderr.decode("UTF-8")
        stderr="\n".join(textwrap.wrap(
                stderr, break_on_hyphens=False, break_long_words=False))

        self.logger.error("ESTSeq.pl failed!")
        self.logger.error("STDERR: {0}".format(stderr))
        self.failed=True
        return

    subprocess.call(['gzip', blatout], shell=False)
    return

#############################################

def doLast(self, masked_fasta=True):
    '''Function to create the .align file. IMPORTANT: the database should be a **collapsed** FASTA file, with only one sequence.
    The name of the two FASTA **must** equate the name of the sequence.'''

    #I have not got around to do this yet.
    raise NotImplementedError

#     if masked_fasta: addition="masked."
#     else: addition=""

#     fasta=os.path.join(dataset, "{0}.{1}fa".format(name, addition))

#     alignOut=os.path.join(dataset, "{0}.{1}.{2}align".format(name, os.path.basename(last_informant), addition)) #The name of the database = name of the species

#     if os.path.exists(alignOut) and os.stat(alignOut).st_size>0:
#         return 0

#     lastOut=os.path.join(dataset, "{0}.{1}.maf".format(name, os.path.basename(last_informant)))
#     if not os.path.exists(lastOut) or os.stat(lastOut).st_size==0:
#         lastOut=open(lastOut, 'w')
#         last_call=subprocess.call(['lastal', last_informant, fasta], shell=False, stdout=lastOut)
#         if last_call!=0:
#             return last_call
#         lastOut.close()
#         lastOut=lastOut.name
        
#     mock_file=os.path.join(dataset, name)
#     ln_call=subprocess.call(['ln', '-s', fasta, mock_file], shell=False)
#     alignOut=open(alignOut, 'w')
#     maf2align_command = subprocess.Popen(['/opt/N-SCAN/bin/maf_to_align.pl',  tempfile.gettempdir(), lastOut, 'A',
#                                          mock_file, last_informant],
#                                         shell=False, stdout=alignOut, stderr=subprocess.PIPE)
#     maf2align_command.communicate()
#     alignOut.close()
#     unlink_call=subprocess.Popen(['unlink', mock_file], shell=False)
#     unlink_call.communicate()
# #    return maf2align_command.returncode
#     return True
    

##############################################
