#!/usr/bin/env python3
#coding: utf_8

import io

class gffLine(object):

    def __init__(self,line, my_line='', header=False):
        '''Object which serializes a GFF line.
		Parameters:
				- _line: the original line
				- _fields: the splitted line
				- chrom: the chromosome
				- source: source field, where the data originated from.
				- feature: mRNA, gene, exon, start/stop_codon, etc.
				- start: start of the feature
				- end: stop of the feature
				- strand: strand of the feature
				- score
				- phase
				- attributes - a dictionary which contains the extra information.
				
				Typical fields in attributes are: ID/Parent, Name'''

        self.attributes=dict()
        self.id=None
        self.parent=None

        if line=='' and my_line!="":
            self._line=my_line
        else:
            self._line=line
        self._fields=line.rstrip().split('\t')
        self.header=header
        # if len(self._fields)!=9:
        #     print(*self._line, file=sys.stderr)

        if self.header or len(self._fields)!=9 or self._line=='':
            self.feature=None
            return

        if len(self._fields)!=9: return None
        self.chrom,self.source,self.feature=self._fields[0:3]
        self.start,self.end=tuple(int(i) for i in self._fields[3:5])

        if self._fields[5]=='.': self.score=None
        else: self.score=float(self._fields[5])

        self.strand=self._fields[6]
        if self.strand=='.': self.strand=None
        assert self.strand in (None,'+','-','?')

        if self._fields[7]=='.': self.phase=None
        else: 
            try: 
                self.phase=int(self._fields[7]); assert self.phase in (0,1,2)
            except: raise

        self._Attr=self._fields[8]

        self.attributeOrder=[]

        for item in [x for x in self._Attr.rstrip().split(';') if x!='']:
            itemized=item.split('=')
            try:
                self.attributes[itemized[0]]=itemized[1]
                self.attributeOrder.append(itemized[0])
            except IndexError:
                pass
#                raise IndexError(item, itemized, self._Attr)
            
        if self.id is None:
            id_key = list(filter(lambda x: x.upper()=="ID", self.attributes.keys()))
            if len(id_key)>0:
                id_key=id_key[0]
                self.id = self.attributes[id_key]
        if self.parent is None:
            parent_key = list(filter(lambda x: x.upper()=="PARENT", self.attributes.keys()))
            if len(parent_key)>0:
                parent_key = parent_key[0]
                self.parent = self.attributes[parent_key]

        assert self.parent is not None or self.id is not None
            
        if "PARENT" in self.attributes and "Parent" not in self.attributes:
            self.attributes['Parent']=self.attributes['PARENT'][:]
            del self.attributes['PARENT']

    @property
    def id(self):
        return self.attributes["ID"]
    
    @id.setter
    def id(self,Id):
        self.attributes["ID"]=Id
        
    @property
    def parent(self):
        return self.attributes["Parent"]
    @parent.setter
    def parent(self,parent):
        self.attributes["Parent"]=parent
    

    def __repr__(self):

        return "\t".join(str(x) for x in [self.chrom, self.source, self.feature, self.start, self.end, self.score, self.strand, self.phase])

    def __str__(self): 
        if not self.feature: return self._line.rstrip()

        if "score" in self.__dict__ and self.score: score=str(self.score)
        else: score="."
        if 'strand' not in self.__dict__ or not self.strand: strand="."
        else: strand=self.strand
        if self.phase!=None: phase=str(self.phase)
        else: phase="."
        attrs=[]
        if self.id is None and self.parent is None:
            raise ValueError("This object has no parent nor ID!\n{0}".format(repr(self)))
        if self.id is not None:
            attrs=["ID={0}".format(self.id)]
            
        if self.parent is not None:
            attrs.append("Parent={0}".format(self.parent))
        for att in self.attributeOrder:
            if att in ["ID","Parent"]: continue
            try: attrs.append("{0}={1}".format(att, self.attributes[att]))
            except: continue #Hack for those times when we modify the attributes at runtime
            
        line='\t'.join(
            [self.chrom, self.source,
             self.feature, str(self.start), str(self.end),
             str(score), strand, phase,
             ";".join(attrs)]
        )
        return line

    def __len__(self):
        if "end" in self.__dict__:
            return self.end-self.start+1
        else: return 0


class GFF3(object):
    def __init__(self,handle):
        if isinstance(handle,io.IOBase):
            self._handle=handle

        else:
            assert isinstance(handle,str)
            try: self._handle=open(handle)
            except: raise ValueError('File not found: {0}'.format(handle))

        self.header=False

    def __iter__(self): return self

    def __next__(self):
        line=self._handle.readline()
        if line=='': raise StopIteration

        if line[0]=="#":
            return gffLine(line, header=True)

        return gffLine(line)
