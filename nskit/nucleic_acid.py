from functools import cached_property
from typing import Optional, Union
from pathlib import Path
import numpy
import numpy as np

from .nucleic_acid_graph import NucleicAcidGraph
from .draw import DrawNA
from .exceptions import InvalidSequence, InvalidAdjacency, InvalidStructure



class NucleicAcid(NucleicAcidGraph, DrawNA):
    __slots__ = ('__name', '__meta')
    
    def __init__(self):
        super().__init__()
        self.__name = None
        self.__meta = None
    
    
    @property
    def name(self) -> str:
        return self.__name or ''
    
    
    @name.setter
    def name(self, name):
        if not isinstance(name, str):
            raise TypeError('NA name must be a string')
        
        name = name.strip()
        if name.startswith('>'):
            raise ValueError("'>' should not be used in name, it'll be automaticaly added during writing")
        self.__name = name
    
    
    @property
    def meta(self) -> dict:
        if self.__meta is None:
            self.__meta = {}
        return self.__meta
        
        
    @property
    def seq(self) -> str:
        return ''.join([nb for _, nb in self._nodes.items()])
    
    
    @cached_property
    def struct(self) -> str:
        return self.assemble_dot_structure()
    
    
    def __str__(self):
        s = f"{self.seq}\n{self.struct}"
        if self.__name:
            s = f"{self.__name}\n{s}"
        return s
    
    
    def __repr__(self):
        return str(self)
    
    
    def __len__(self):
        return len(self._nodes)
    
    
    def __eq__(self, other):
        if isinstance(other, str):
            return self.seq == other.strip().upper()
        
        elif isinstance(other, NucleicAcid):
            return self.seq == other.seq
        
        else:
            raise TypeError(f"Can not compare NucleicAcid and {type(other)}")
    
    
    @classmethod
    def from_adjacency(cls, adj: numpy.array, /, *, 
                       seq: Optional[str] = None,
                       name: Optional[str] = None,
                       meta: Optional[dict] = None,
                       upper_sequence: bool = True, 
                       fix_sharp_helixes: bool = False, 
                       trust_adj: bool = False
                      ) -> 'NucleicAcid':
        """
        Create NucleicAcid from adjacency matrix.

        :param adj: numpy array adjacency of 0 and 1.
        :param seq: na sequence.
        :param name: na name.
        :param meta: dictionary of meta information convertable to string.
        :param upper_sequence: upper sequence characters. Default - True.
        :param fix_sharp_helixes: remove bond between neighboring nbs to fix sharp helixes. Default - False.
        :param trust_adj: whether to skip adjacency validation. Validation has O(N^2) time complexity. Default - False.

        :return: NucleicAcid object.
        """ 
        
        if len(adj.shape)!=2 or adj.shape[0]!=adj.shape[1]:
            raise InvalidAdjacency(f"Adjacency must be a square matrix, got shape: {adj.shape}")
            
        if seq is None: 
            seq = 'N'*adj.shape[0]
        else:
            if not seq.isalpha():
                raise InvalidSequence(f"Sequence must contain only alphabetic characters")
            
            if upper_sequence: 
                seq = seq.upper()

            if len(seq)!=adj.shape[0]:
                raise InvalidAdjacency(f"Adjacency shape and sequence length must be equal, got adj: {adj.shape[0]} and seq: {len(seq)}")
        
        if not trust_adj:
            if np.sum((adj!=0)*(adj!=1)) > 0:
                raise InvalidAdjacency(f"Adjacency must contain only 0 and 1")
            
            if not np.array_equal(adj, adj.T):
                raise InvalidAdjacency(f"Adjacency must be symmetric")
            
            if np.sum(adj.sum(axis=-1) > 1) > 0:
                raise InvalidAdjacency(f"Several complementary bonds for one nucleic base is ambiguous")
        
        # create graph
        na = cls()
        if name: na.name = name
        if meta: na.meta.update(meta)
            
        for i, nb in enumerate(seq):
            _ = na.add_node(nb)
            if i>0:
                na.add_bond(i-1, i, 0)
            
        adj = np.triu(adj, 1) # mask diagonal and lower triangle
        vec = np.argmax(adj, axis=-1)
        for o, e in enumerate(vec):
            if e==0: 
                continue

            if abs(o-e)==1:
                if fix_sharp_helixes:
                    continue
                raise InvalidStructure(f"Sharp helix between nbs {o} and {e}, use fix_sharp_helixes=True to omit such bonds")
                
            na.add_bond(o, e, 1)
            
        return na
    

    def to_dot(self, path: Union[str, Path], *, 
               append: bool = False,
               write_struct: bool = True, 
               write_meta: bool = True
               ):
        
        from .io import DotWrite
        
        with DotWrite(path, append=append) as w:
            w.write(self, write_struct=write_struct, write_meta=write_meta)

    
    def to_fasta(self, path: Union[str, Path], *, 
                 append: bool = False, 
                 ):
        from .io import FastaWrite
        
        with FastaWrite(path, append=append) as w:
            w.write(self)


    def to_bpseq(self, path: Union[str, Path], *, 
                 write_meta: bool = True
                 ):
        from .io import bpseqWrite
        
        with bpseqWrite(path) as w:
            w.write(self, write_meta=write_meta)