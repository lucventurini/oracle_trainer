#!/usr/bin/env python

from Oracle_utils.interface import Interface
import sys,subprocess,os,logging
import re, tempfile
from myRecords import GTF
from Bio import SeqIO
import textwrap


class Initializer(Interface):

    '''This class prepares the basic files for each gene model:
    - fasta files (masked and unmasked)
    - gff3 files
    - gtf files
    It also performs various quality checks on the GTF to ensure everything is set up properly.'''

    def __init__(self, interface=None):
        
        '''Initializer function.'''


        if interface is None:
            raise ValueError("No interface provided!")

        self.inherit(interface) #Inherit logger, datafolder, etc.
        self.define_step()


        ##Additional parameters
        self.chrom=None
        self.start=float("Inf")
        self.end=float("-Inf")
        self.strand=None
        self.failed=False

#############################

    def define_step(self):
            self.step="INITIALIZATION"

############################

    def __call__(self):

        '''This function calls the different methods to create the necessary files.
        If anything goes amiss, it returns immediately and reports the failure.
        Resolution order:
        - checkProtein
        - createGtf
        - createSequences
        - validateGtf
        - gtf2gff3'''

        order=["checkProtein", "createGtf", "createSequences",
               "validateGtf", "gtf2gff3"]

        for method in order:
            getattr(self, method)()
            if self.failed==True:
                break

        super(self.__class__, self).__call__()
            
##############################################

    def checkProtein(self):
        '''Simple sequence checker. It ensures that no masked bases are present and that the first AA is a methionine.'''

        self.logger.info("Checking FASTA file.")
        if self.record.seq[0]!="M" or "X" in str(self.record.seq):
            self.logger.critical("Irregular protein sequence. Exiting.")
            self.failed=True


#############################################

    def createFasta(self, masked=False):

        '''This function uses SeqIO.index_db to recover the slice pertaining to the gene model.'''

        self.logger.debug("Recovering the sequence from chromosome {0} from {1} to {2}".format(self.chrom, self.start, self.end))

        if masked:
            out = os.path.join(
                self.dataFolder, "{0}.masked.fa".format(self.new_name))
            index=SeqIO.index_db(self.masked_fasta)
            self.masked_fasta = out

        else:
            out=os.path.join(
                self.dataFolder, "{0}.fa".format(self.new_name))
            index=SeqIO.index_db(self.genomic_fasta)
            self.fasta = out

        if self.chrom not in index:
            self.failed=True
            self.logger.critical("Chromosome {0} not in index".format(self.chrom))
            return

        length=len(index[self.chrom])
        #Define boundaries
        start=max(0, self.start-1-self.flank)
        end=min(length-1, self.end-1+self.flank)
        
        if start>=end:
            self.logger.critical("Invalid values for this chromosome: {start}-{end}".format(start, end))
            self.failed=True
            return

        seq=index[self.chrom][start:end]
        with open(out,'w') as out:
            seq.description = ''
            seq.id = self.new_name
            print(seq.format('fasta'), file=out, end='')


    def createSequences(self):
        '''Thin wrapper for createFasta.'''

        for flag in [False,True]:
            self.createFasta(masked=flag)
            if self.failed: break

        return

#######################################################################
    def createGtf(self):

        '''This function creates the GTFs (with an without UTRs).
        Contestually, it will also perform sanity checks such as verifying that the model has start and stop codons.'''

        self.logger.info("Beginning the creation of the GTFs.")

        #Set the output filenames
        self.gtf=os.path.join(
            self.dataFolder, "{0}.gtf".format(self.new_name))
        self.gtf_no_utr = os.path.join(
            self.dataFolder, "{0}.noutr.gtf".format(self.new_name))

        if os.path.exists(self.gtf) and os.path.getsize(self.gtf)>0:
            if os.path.exists(self.gtf_no_utr) and os.path.getsize(self.gtf_no_utr)>0:
                self.logger.warn("GTFs already present. Skipping.")
                return

        gtf_lines = []

        self.logger.debug("Parsing the original GTF.")
        for gtfline in GTF.GTF(open(self.gengtf)):
            if gtfline.transcript==self.name:
                if not self.chrom:
                    self.chrom=gtfline.chrom
                self.start=min(self.start, gtfline.start)
                self.end=max(self.end, gtfline.stop)
                if not self.strand: self.strand=gtfline.strand
                gtf_lines.append(gtfline)

        #Sanity checks
        if not self.chrom:
            self.logger.error("No chromosome defined!")
            self.failed=True
        if not len([x for x in gtf_lines if x.feature=="start_codon"])==1:
            self.logger.error("No start codon defined!")
            self.failed=True
        if not len([x for x in gtf_lines if x.feature=="stop_codon"])==1:
            self.logger.error("No stop codon defined!")
            self.failed=True
        positions = [ (line.start, line.stop) for line in 
                       filter(lambda gtf_line: gtf_line.feature=="CDS", gtf_lines)]

        if len(positions)!=len(set(positions)):
            self.logger.error("Duplicated exons in record!")
            self.failed=True

        #If we have failed at some point, exit
        if self.failed==True:
            return

        self.offset = max(0, self.start - self.flank-1)

        gtf=open(self.gtf,'w')
        gtf_no_utr=open(self.gtf_no_utr,'w')
        self.logger.debug("GTF no UTR: {0}, handle {1}".format(self.gtf_no_utr, gtf_no_utr))

        for line in gtf_lines:
            line.chrom = self.new_name
            line.start -= self.offset
            line.stop -= self.offset
            if line.feature in ("CDS", "start_codon", "stop_codon"):
                print(line, file=gtf_no_utr)
            print(line, file=gtf)

        gtf.close()
        gtf_no_utr.close()
        self.logger.info("Finished preparing the GTFs.")
        return
            
#############################################################

    def validateGtf(self):

        '''This method uses the EVAL_GTF script validate_gtf.pl to check the GTF for consistency.'''

        validation_file = open(os.path.join(
                self.dataFolder, "validation.txt"), 'w')


        self.logger.info("Checking the GTF for consistency.")
        
        command_line=[os.path.join(os.environ['EVAL_GTF'],  'validate_gtf.pl')]
        command_line+=[self.gtf,self.fasta]

        validate_gtf = subprocess.Popen( command_line,
                                         shell=False, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
        
        stdout, stderr = validate_gtf.communicate()
        if validate_gtf.returncode not in (None,0):
            self.logger.critical("Validate_gtf failed. STDERR:")
            stderr="\n".join(textwrap.wrap(stderr.decode("UTF-8")))


        stdout = stdout.decode("UTF-8")
        lines=stdout.split("\n")

        warnings = list(filter( lambda line: re.search("Warnings", line), lines))
        if len(warnings)>0:
            if self.only_canonical:
                self.failed=True
                return
            else:
                lines=list(filter(lambda line: line!='', lines))
                war_position = lines.index("Warnings encountered:")+1
                stat_pos = lines.index("Statistics:")
                warning_lines = lines[ war_position+1:stat_pos ]
                
                self.logger.debug("Warning lines:"+"\n".join(warning_lines))

                if len(list(filter(lambda line: "non-canonical splice" not in line, warning_lines)))>0:
                    self.failed=True
                    return
                    
        print(stdout, file=validation_file)
        return

###########################

    def gtf2gff3(self):
        '''This method makes use of the gtf2gff3 script from GAL,
        provided by the Sequence Ontology Consortium (SO).
        Webpage: http://www.sequenceontology.org/software/GAL.html'''

        self.gff3 = os.path.join( self.dataFolder, "{0}.gff3".format(self.new_name))

        if os.path.exists( self.gff3 ) and os.path.getsize(self.gff3)>0:
            self.logger.warn("GFF3 already present. Skipping.")
            return

        out=open(self.gff3, 'w')
        gff3conversion = subprocess.Popen( ['gtf2gff3', self.gtf],
                                           shell=False,
                                           stdout=out,
                                           stderr=subprocess.PIPE)

        stdout,stderr = gff3conversion.communicate()
        if gff3conversion.returncode not in (0, None):
            stderr=stderr.decode("UTF-8")
            os.path.remove(self.gff3)
            self.logger.critical("Conversion to GFF3 failed! STDerr:")
            self.logger.critical("\n".join( textwrap.wrap(
                        stderr, break_on_hyphens=False, break_long_words=False)))
            self.failed=True
            return

############################
