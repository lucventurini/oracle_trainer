#!/usr/bin/env python3


import sys,re, itertools,os,argparse
from operator import itemgetter,attrgetter,add
from Bio import SeqIO, Seq,SeqRecord
from math import log,e
from scipy import mean,std,array
import argparse
import scipy as sp
from functools import reduce
sp.seterr('raise')


'''This function is a port of the MarkovMatrices*awk functions. It takes as input the order of the Markov model AND the fasta file, and it constructs the relative prob. matrix. These are later used to construct the true MM with the stubs pro2log*awk'''

class Markov(object):

    def __init__(self, observed=None, background=None,
                 transition=None, initial=None,
                 signal_probs=None,
                 retrieve=False, signal=None, length=0,
                 signal_maximum=1, cds_maximum=5, window_size=5, min_zscore = 0.2,
                 cutoff=-5, folder="."):
        '''Init function. It takes as input either the background + observed FASTA files or, in alternative, the pre-processed counts.
        Signal must be in (None,"Acceptor", "Donor", "Start", "Stop"); if None, it is assumed we are analyzing CDS vs. intron.
        Arguments:
            observed: EITHER the indexed FASTA of observed signals, OR the observed frequencies and ratios for signal regions.
            background: the indexed FASTA of background information.
            transition: the observed transition frequences for CDS kmers.
            raw_counts: whether we are retrieving raw counts from the files. Not sure I have implemented it.
            signal: {None (cds and introns), Start, Donor, Acceptor}
            length: length of the sequences (for signals)
            top_one: for Acceptor and Donor signals, whether the maximum order of the matrix must be one (True) or two (False).
            min_zscore: minimum averaged Z-score for a position to be deemed significant.
                        Default: 0.2 (as in the publication for SPARCS, see doi:10.1093/nar/gkt461 )
            window: Lenght of the overlapping windows for the calculation of the averaged Z-Scores. Default: 5.
        '''
        #Initialize the count dictionary
        self.cutoff=cutoff
        self.counts={}
        self.probabilities={}
        self.counts['background']={}
        self.counts['observed']={}
        self.signals={"Acceptor": "AG",
                      "Donor": "GT",
                      "Start": "ATG",}

        #Store various parameters
        assert signal in (None,"Acceptor", "Donor", "Start")
        self.signal=signal
        self.length=length
        self.window_size = window_size
        self.min_zscore = min_zscore
        self.signal_maximum = signal_maximum
        self.folder=folder
        self.cds_maximum=cds_maximum

        ##Files
        self.observed=observed
        self.background=background
        self.initial=initial
        self.transition=transition
        self.signal_probs=signal_probs


        if retrieve:
            assert (self.signal_probs!=None and self.signal!=None) or (self.initial!=None and self.transition!=None), \
                (self.signal_probs, self.signal, self.initial, self.transition)
            self.retrieveFrequencies()
        else:
            assert self.background!=None, "A FASTA file for both observed and background sequences must be provided!"
            self.calculateOrder(background, observed, cds_maximum=cds_maximum, signal_maximum=signal_maximum)
            if self.signal:
                try:
                    if self.length==0: self.length=len(list(background.keys())[0])
                except IndexError:
                    raise IndexError("Index out of range", background, list(background.keys()))
                self.counts['background']=self.signal_counts(background)
                self.counts['observed']=self.signal_counts(observed)
            else:
                '''This part of the function treats the CDS/Intron case.'''
                self.counts['background']['initial']={}
                self.counts['background']['transition']={}
                self.counts['background']['initial'],self.counts['background']['transition'],self.background_total= self.cds_counts(background, intron=True)
                
                self.counts['observed']['initial']={}
                self.counts['observed']['transition']={}
                self.counts['observed']['initial'],self.counts['observed']['transition'],self.observed_total= self.cds_counts(observed, intron=False)

            #Remove the dictionaries
            del observed
            del background
            
            #Calculate the probabilities
            self.compute()

#########################

    def print_markov(self):
        '''Mock method which will call the right formatting functions.'''

        if self.signal:
            for line in self.printSignal(): yield line
        else:
            for line in self.printCds(): yield line


    def printSignal(self):
        ''''This method elaborates and print the signal 
        output matrix for the geneid file.'''
#        self.signalBoundaries()
#        self.compute()
        self.signalBoundaries()
#        print(self.matrixStart, self.matrixStop, file=sys.stderr)

        keys=sorted(sorted(self.probabilities, key=itemgetter(1)), key=itemgetter(0)) 
        max_profile_length=max([ key[0]-self.matrixStart+1 for key in keys])

        profileLength=min(self.matrixStop+1-self.matrixStart, max_profile_length)

    
        yield "{0}_profile".format(self.signal)
        yield " ".join([str(piece) for piece in [profileLength, self.offset, self.cutoff, self.order]])


        for key in keys:
            if key[0] not in list(range(self.matrixStart, self.matrixStop+1)): continue
            pos=key[0]-self.matrixStart+1
            kmer=key[1]
            probability=self.probabilities[key][-1]

            line=[str(k) for k in 
                  [pos, kmer, probability]
                  ]
            yield " ".join([str(k) for k in line])

        return

    def printCds(self):
        '''This function prints the Initial and Transition matrices into the
        geneid parameter file.'''

#        self.compute()
        yield str(self.order)
        initial_count=-1
        yield "Markov_Initial_probability_matrix"
        
        #Keys: frame, kmer
        keys=sorted(sorted(self.probabilities['initial'], key=itemgetter(0)), key=itemgetter(1))

        for key in keys:
            if key[0]==0: initial_count+=1
            line=[key[1], initial_count, key[0], self.probabilities['initial'][key]]
            yield " ".join([str(k) for k in line])

        yield ""

        transition_count=-1
        yield "Markov_Transition_probability_matrix"
        keys=sorted(sorted(list(self.probabilities['transition'].keys()), key=itemgetter(0)),key=itemgetter(2))
        for key in keys:
            if key[0]==0: transition_count+=1
            line=[key[2], transition_count, key[0], self.probabilities['transition'][key]]
            yield " ".join([str(k) for k in line])
            

#########################

    def compute(self):
        '''Mock method to access the right computation method.'''

        if self.signal:
            self.probabilities=self.compute_signal()
        else:
            self.frequency={}
            self.frequency['initial']={}
            self.frequency['transition']={}
            assert isinstance(self.background_total, int), "Compute wrong"
            assert isinstance(self.observed_total, int), "Compute wrong"
            self.frequency['initial']['background'],self.frequency['transition']['background'] = self.cds_frequency('background', self.background_total, intron=True)
            self.frequency['initial']['observed'], self.frequency['transition']['observed'] = self.cds_frequency('observed', self.observed_total)
            self.probabilities['initial'],self.probabilities['transition']=self.compute_cds()


#########################

    def cds_counts(self, fasta, intron=False):
        '''This method reimplements the first part of the MarkoMatrices.awk script.
        In particular, it calculates the counts of each possible kmer 
        (derived by the order of the matrix) in the various input sequences.'''
       
        pcount=1/4**(self.order+1) #Pseudocount for the matrices
        initial={}
        transition={}

        nucs=["A","C","G","T"]
        ini_kmers=["".join(el) for el in sorted(set(itertools.combinations(nucs*self.order,self.order)))]
        trans_kmers=["".join(el) for el in sorted(set(itertools.combinations(nucs*(self.order+1),self.order+1)))]
        if not intron:
            for iii in range(3):
                for kmer in ini_kmers: initial[(iii,kmer)]=pcount*4
                for kmer in trans_kmers: transition[(iii,kmer)]=pcount
        else:
            for kmer in ini_kmers: initial[kmer]=pcount*4
            for kmer in trans_kmers: transition[kmer]=pcount

        ##Initialization done
        total_length=0
        for sequence in fasta:
            num=fasta[sequence]
            total_length+=(len(sequence)-self.order)*num #This will be useful at the end
            if not intron:

                frames = [iii%3 for iii in range(len(sequence)-self.order) ]
                ini_kmers = zip(
                    frames,
                    [reduce(add, sequence[pos:pos+self.order]) for pos in range(len(sequence)-self.order)]
                    )

                trans_kmers = zip(
                    frames,
                    [reduce(add, sequence[pos:pos+self.order+1]) for pos in range(len(sequence)-self.order-1)]
                    )
            else:
                #For introns we do NOT consider the frame!

                ini_kmers = [reduce(add, sequence[pos:pos+self.order]) for pos in range(len(sequence)-self.order)]
                trans_kmers = [reduce(add, sequence[pos:pos+self.order+1]) for pos in range(len(sequence)-self.order-1)]

            for key in ini_kmers:
                self.add_num(initial, key, num)
            for key in trans_kmers:
                self.add_num(transition, key, num)

            # for iii in xrange(len(seq)-self.order):
            #     ini_kmer=seq[iii:iii+self.order]
            #     trans_kmer=seq[iii:iii+self.order+1]
            #     if re.search("[^ACGT]", ini_kmer) or re.search("[^ACGT]",trans_kmer): continue

            #     if not intron:
            #         frame=iii%3
            #         initial[(frame, ini_kmer)]+=1
            #         transition[(frame, trans_kmer)]+=1
            #     else:
            #         initial[ini_kmer]+=1
            #         transition[trans_kmer]+=1

        return initial, transition, total_length

###########################

    def cds_frequency(self, key, total_length, intron=False):
        ''''Second part of the MarkovMatrices script.
        It is used to calculate the probabilities of each kmer
        and the probability of each possible transition,
        according to the data.
        If we are working with CDSs, it also indexes the kmers by their frame.
        The latter aspect is controlled by the keyword "intron".
        Arguments:
            - key: the key to look for in dictionary ("background" or "observed")
            - total_length: the amount of sequence in the data. Necessary to calculate the frequencies!
            - intron [keyword]: whether we have to index also by frame, or not.'''

        initial=self.counts[key]['initial']
        transition=self.counts[key]['transition']
        
        if intron: factor=1
        else: factor=3
        assert isinstance(total_length, int), "Wrong frequency"
        total_length+=factor*(4**self.order)

        initial_probs=dict().fromkeys(initial)
        for key in initial:
            initial_probs[key]=initial[key]/(total_length/factor)

        transition_probs={}
        for key in transition:
            if not intron:
                frame=key[0]
                kmer=key[1]
                ini_kmer=key[1][:-1]
                transition_probs[(frame,
                                  ini_kmer,
                                  kmer)]=transition[key]/initial[(frame,ini_kmer)]
            else:
                kmer=key
                ini_kmer=key[:-1]
                transition_probs[(ini_kmer, kmer)]= \
                    transition[kmer]/initial[ini_kmer]

        return initial_probs, transition_probs

##############################

    def add_num(self, dictionary, key, num):
        '''Quick function to add any number to the kmer occurrence'''
        try: dictionary[key]+=num
        except KeyError:
            pass

    def signal_counts(self, fasta):
        '''This method copies the frequency.awk script.
        It is used to count the occurrences of the kmers inside the signal frequencies.'''
        #Initialize the dictionary of the occurrences

        occurrences=dict()
        nucs=["A","C","G","T"]
        k=self.order+1
        kmers=["".join(el) for el in sorted(set(itertools.combinations(nucs*k,k)))]
        count=0
        for iii in range(1,self.length+1-self.order,1):
            for kmer in kmers:
                occurrences[(iii,kmer)]=0
        
        for record in fasta:
            num = fasta[record]
            
            #Use reduce to create the kmer-list and izip to link it to the position
#            print(record, file=sys.stderr)
            kmers=zip( 
                range(1,len(record)-self.order+1), 
                [reduce(add, record[pos:pos+self.order+1]) for pos in range(len(record)-self.order)],
                )

            array_length=len(record)-self.order
            list(map(self.add_num, (occurrences,)*array_length, kmers, (num,)*array_length))
            count+=num
                        
            # for l in xrange(1,self.length+1-self.order,1):
            #     kmer=str(seq)[l-1:l-1+k].upper()
            #     if len(kmer)!=self.order+1: continue #Bug
            #     if re.search("[^ACGT]", kmer): continue
            #     try: occurrences[(l,kmer)]+=1
            #     except:
            #         raise ValueError, (record, "\n", l, kmer)
        
        for key in occurrences:
            try:
                occurrences[key]=(occurrences[key], occurrences[key]/count)
            except ZeroDivisionError:
                occurrences[key]=(occurrences[key], 0)
        return occurrences


#############################

    def compute_signal(self):

        '''This method re-implements the aux.awk script, and 
        outputs a dictionary with the probabilities.'''

        center=int(self.length/2)+1
        observed=self.counts['observed']
        background=self.counts['background']
        
#        sigRange=range(center, center+1+int(len(self.signals[self.signal]))) #These are the conserved positions
        if self.signal=="Start":
            #For the start signal, the 0s must stop right after we have determined the ATG
            sigRange=list(range(
                center-self.order,
                center))
        else:
            sigRange=list(range(
                center-self.order,
                center+len(self.signals[self.signal])-1+self.order-1))


        probabilities={}

        table_out=open(os.path.join(self.folder, "{name}_freqs.txt".format(name=self.signal.lower())), 'w')

        for pos in range(1,self.length+1):
            for key in [x for x in list(self.counts['observed'].keys()) if x[0]==pos]:
                val=[observed[key][1], background[key][1]]
                
                #Inside signal range
                if key[0] in sigRange:
                    #Signal
                    if observed[key][1]!=0:
                        val.append(0)
                    #Non-signal
                    else: val.append(-9999)

                #Outside of the signal range
                elif observed[key][1]!=0 and background[key][1]!=0:
                    val.append( log(observed[key][1]/background[key][1],e) )
                #Outside of the signal range, no data
                else: val.append(-9999)
                probabilities[key]=val

                table_line=list(key)+val
                print(*table_line, sep="\t", file=table_out)

        return probabilities

############################

    def compute_cds(self):
        '''This method reimplements the pro2log scripts.
        Its purpose is to calculate the logarithm of the ratio of 
        observed vs. background k-mer counts in the initial
        and transition matrices (i.e. ratio of the frequency
        of each kmer in observed vs. background populations
        and frequency of each possible transition.'''

        initial_out=open(os.path.join(self.folder, "cds_initial_freqs.txt"), 'w')
        initial_probabilities={}
        for key in self.frequency['initial']['observed']:
            ##Structure of the keys: (frame, kmer)
            frame, kmer= key
            obs=self.frequency['initial']['observed'][key]
            back=self.frequency['initial']['background'][kmer]
            initial_probabilities[key] = log(obs/back, e)
            init_line=list(key)+[obs,back,initial_probabilities[key]]
            print(*init_line, sep="\t", file=initial_out)

        initial_out.close()

        transition_probabilities={}
        transition_out=open(os.path.join(self.folder, "cds_transition_freqs.txt"), 'w')
        for key in self.frequency['transition']['observed']:
            #Structure of the keys: (frame, ini_kmer, trans_kmer)
            frame, ini_kmer, trans_kmer = key
            obs=self.frequency['transition']['observed'][key]
            back=self.frequency['transition']['background'][(ini_kmer, trans_kmer)]
            transition_probabilities[key]=log(obs/back, e)
            trans_line=list(key)+[obs, back, transition_probabilities[key]]
            print(*trans_line, sep="\t", file=transition_out)

        transition_out.close()

        return initial_probabilities, transition_probabilities

############################

    def signalBoundaries(self):
        '''This function calculates the boundary of the region with significant information around the signal considered.
        It performs the calculation by computing the Kullback-Leibler entropy for each position and subsequently smoothing
        with overlapping windows of a given length (parameter window_size.)
        The significant positions are those whose average across the overlapping windows is greater than a set value.
        This method has been inspired by the similar algorithm in SPARCS ( doi:10.1093/nar/gkt461 )'''

        assert self.signal!=None
        
        center=int(self.length/2)+1
        sigRange=list(range(center, center+len(self.signals[self.signal])+1))
        
        if self.signal=="Start":
            sigStart=sigRange[0]
            minimum=sigStart-2
            matrixStart=sigRange[0]-2
            sigStop=sigStart+2
        elif self.signal=="Donor":
            sigStart=sigRange[0]-2
            minimum=sigStart-1
            matrixStart=sigRange[0]-3
            sigStop=sigStart+1
        elif self.signal=="Acceptor":
            sigStart=sigRange[-1]
            matrixStart=sigRange[0]-1
            minimum=matrixStart-1
            sigStop=sigStart+1

        #Calculate entropies for each position
        entropies=dict()

        for pos in range(1,self.length+1):
            prob_keys=[prob_key for prob_key in self.probabilities if prob_key[0]==pos]
            vals=[]
            for key in prob_keys:
                #Kullback-Leibler Entropy is calculated as SUM(obs_freq * ln(obs_freq/back_freq) ) 
                vals.append(self.probabilities[key][0]*self.probabilities[key][-1])
            entropies[pos]=sum(vals)

        mean_entropy=mean(list(entropies.values()))
        std_entropy=std(list(entropies.values()))

        #Create the windows for the smoothing

        windows=dict()
        for num in range(1,self.length+1 - self.window_size):
            key=tuple(range(num, num+1+self.window_size))
            vals=[]
            for pos in key:
                vals.append(entropies[pos])
            
            #Zcore = (Mean(sample) - Mean(population))/(STD(population))
            zscore = ( mean(vals)-mean_entropy )/ std_entropy #IMPLEMENT
            windows[key]=zscore

        smoothed_values=dict()
        for num in range(1,self.length+1):
            my_windows = [key for key in list(windows.keys()) if num in key]
            smoothed_values[num]=mean([windows[key] for key in my_windows])


        values_to_retain=[num for num in list(smoothed_values.keys()) if smoothed_values[num] >= self.min_zscore]


        #The stop and start of the signal MUST be in the values
        for pos in range(sigStart,sigStop+1):
            if pos not in values_to_retain:
                values_to_retain.append(pos)

        values_to_retain=sorted(values_to_retain)


        entropy_out = open(os.path.join(self.folder, "{name}_entropy.txt".format(name=self.signal.lower())), 'w')
        for pos in sorted(smoothed_values.keys()):
            print(pos, smoothed_values[pos], sep="\t", file=entropy_out)
        entropy_out.close()


        #Find the interval
        
        sigStartPos = values_to_retain.index(sigStart)
        sigStopPos = values_to_retain.index(sigStop)
        matrixStartPos=sigStartPos
        matrixStopPos=sigStopPos

        for num in range(sigStartPos-1,0,-1):
            my_pos = values_to_retain[num]
            previous = values_to_retain [matrixStartPos]
            if my_pos+1==previous:
                matrixStartPos=num
            else:
                break

        for num in range(sigStopPos+1,len(values_to_retain)):
            my_pos = values_to_retain[num]
            previous = values_to_retain[matrixStopPos]
            if my_pos == previous+1:
                matrixStopPos = num
            else:
                break

        matrixStart=min(minimum, values_to_retain[matrixStartPos])
        matrixStop=max(sigStop, values_to_retain[matrixStopPos] )

        self.matrixStart=matrixStart
        self.matrixStop=matrixStop
        if self.signal=="Start":
            self.offset=center-matrixStart #Offset=last non-CDS position
        elif self.signal=="Acceptor":
            self.offset=center+2-matrixStart #Offset=last non-CDS position
        elif self.signal=="Donor":
            self.offset=center-1-matrixStart #Offset=second-to-last CDS position
        self.offset-=self.order
        return 


#################################    

    def retrieveFrequencies(self):
        ''''Mock method to call the right retrieval function'''

        if not self.signal:
            self.retrieveCdsFrequencies()
        else:
            self.retrieveSignalFrequencies()

    def retrieveCdsFrequencies(self):
        assert self.initial!=None and  os.path.exists(self.initial)
        assert self.transition!=None and  os.path.exists(self.transition)

        self.probabilities['initial']=dict()
        self.probabilities['transition']=dict()

        ini_length = 0

        for line in open(self.initial):
            frame,kmer,obs,background,ratio = line.rstrip().split()
            frame=int(frame)
            if len(kmer)>self.cds_maximum:
                raise ValueError("Pre-calculated CDS kmers are greater than the requested maximum!")
            ini_length = len(kmer)
            self.order=len(kmer)-1
            obs,background,ratio=float(obs),float(background),float(ratio)
            #Here, we keep only the ratio - we are not interested in calculating the entropy
            self.probabilities['initial'][tuple( [frame, kmer] ) ] = ratio

        for line in open(self.transition):
            frame, ini_kmer, trans_kmer, obs, background, ratio = line.rstrip().split()
            frame=int(frame)
            if len(ini_kmer)>self.cds_maximum:
                raise ValueError("Pre-calculated CDS kmers are greater than the requested maximum!")
            if len(ini_kmer)+1!=len(trans_kmer):
                raise ValueError("Error in the transition table!")
            key=tuple( [frame, ini_kmer, trans_kmer] )
            #Here, we keep only the ratio - we are not interested in calculating the entropy
            self.probabilities['transition'][key] = ratio

        self.order=ini_length


    def retrieveSignalFrequencies(self):
        assert self.signal_probs!=None and os.path.exists(self.signal_probs)

        kmer_length=0
        for line in open(self.signal_probs):
            pos,kmer,obs,background,ratio = line.rstrip().split()
            pos=int(pos)
            if len(kmer)>self.signal_maximum+1:
                raise ValueError("Pre-calculated signals are greater than the requested maximum!")
            kmer_length=len(kmer)
            obs,background,ratio=float(obs),float(background),float(ratio)
            self.probabilities[ tuple( [pos, kmer] ) ] = [ obs, background, ratio ]
        self.order=kmer_length-1
        self.length=max( [ key[0] for key in self.probabilities] )+self.order


#################################

    ####These are generic utilities

    def calculateOrder(self, background, observed, cds_maximum=5, signal_maximum=1):
        '''This method calculates the maximum order of the computable Markov model.
        Source: http://genome.crg.es/software/geneid/docs/training/html/gif.html'''
        
        ## I really, really, REALLY would like to find a better source for this method.
        self.order=0
        #Now I have a dictionary of the form { (seq, num_seq) }
        backTotal=sum([len(record)*background[record] for record in background])
        observedTotal=sum([len(record)*observed[record] for record in observed])
        
        while backTotal > 30*4**(self.order+1) and observedTotal > 90*4**(self.order+1):
            self.order+=1

        self.order=max(1,self.order-1)
        if self.signal:
            self.order=min(signal_maximum,self.order)
        else:
            self.order=min(cds_maximum, self.order)
        
        return

#################################

    def expr_format(self, signal, n):
        '''This function prepares the regular expression string needed for pattern matching.'''

        before = max(0, -n)
        sigStart = max(0, n)
        sigEnd = min(len(signal), self.order+1+n)
        after = max(0, self.order+1 - len(signal)+n)
        regex="[ACGT]{{{0}}}{1}[ACGT]{{{2}}}".format(before, 
                                                       signal[sigStart:sigEnd], 
                                                       after )

        return  regex

def main():

    parser=argparse.ArgumentParser("Test the module.")
    parser.add_argument('--observed', type=argparse.FileType('r'), default=None)
    parser.add_argument('--background', type=argparse.FileType('r'), default=None)
    parser.add_argument('--initial',  type=argparse.FileType('r'), default=None)
    parser.add_argument('--transition',  type=argparse.FileType('r'), default=None)
    parser.add_argument('--signal_probs',  type=argparse.FileType('r'), default=None)
    parser.add_argument('--signal', default=None, choices=["Acceptor", "Donor", "Start"])
    parser.add_argument('--min_zscore', default=0.2, type=float)
    parser.add_argument('--window_size', default=5, type=int)
    parser.add_argument('--cds_maximum', type=int, default=5)
    args=parser.parse_args()

    if args.observed!=None:
        assert args.background!=None
        args.observed.close()
        args.observed=SeqIO.index(args.observed.name,'fasta')
        args.background.close()
        args.background=SeqIO.index(args.background.name,'fasta')
        markov=Markov(observed=args.observed, background=args.background, signal=args.signal, cds_maximum=args.cds_maximum)

    elif args.signal_probs:
        args.signal_probs.close()
        args.signal_probs=args.signal_probs.name
        markov=Markov(signal_probs=args.signal_probs, retrieve=True, signal_maximum=2, signal=args.signal, 
                      min_zscore = args.min_zscore, window_size=args.window_size)

    elif args.initial:
        args.initial.close()
        args.initial=args.initial.name
        assert args.transition!=None
        args.transition.close()
        args.transition=args.transition.name
#        print(args.cds_maximum, file=sys.stderr)
        markov=Markov(transition=args.transition, initial=args.initial, retrieve=True, cds_maximum=args.cds_maximum)
    else:
        parser.print_usage()
        sys.exit(0)

    for line in markov.print_markov(): print(line)

if __name__=='__main__': main()
