from functools import cached_property
from typing import Optional, Union
from pathlib import Path
import numpy
import numpy as np

from .nucleic_acid_graph import NucleicAcidGraph
from ..draw import DrawNA
from ..exceptions import InvalidSequence, InvalidAdjacency, InvalidStructure



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
        return ''.join(self._nodes)
    
    
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
        """
        Two NucleicAcids are equal if nb sequences and all complementary bonds are equal.
        """
        if not isinstance(other, NucleicAcid):
            raise TypeError(f"NucleicAcid can be compared only with another NucleicAcid, got {type(other)}")
        
        if len(self)!=len(other):
            return False
        
        if self.seq!=other.seq:
            return False
        
        if len(self.pairs)!=len(other.pairs):
            return False
        
        return all([p1==p2 for p1, p2 in zip(self.pairs, other.pairs)])
    

    def __hash__(self) -> int:
        """
        NucleicAcid's hash is computed only by sequence and complementary bonds.
        """
        return hash((self.seq, self.pairs))
    
    
    @classmethod
    def from_adjacency(cls, adj: numpy.array, /, *, 
                       seq: Optional[str] = None,
                       name: Optional[str] = None,
                       meta: Optional[dict] = None,
                       upper_sequence: bool = False,
                       trust_adj: bool = False
                      ) -> 'NucleicAcid':
        """
        Create NucleicAcid from adjacency matrix.

        :param adj: numpy array adjacency of 0 and 1.
        :param seq: na sequence.
        :param name: na name.
        :param meta: dictionary of meta information convertable to string.
        :param upper_sequence: upper sequence characters. Default - False.
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
            _ = na._add_node(nb)
            
        adj = np.triu(adj, 1) # mask diagonal and lower triangle
        vec = np.argmax(adj, axis=-1)
        for o, e in enumerate(vec):
            if e==0: 
                continue
                
            na._add_bond(o, e)
            
        return na
    

    def to_dot(self, path: Union[str, Path], *, 
               append: bool = False,
               write_struct: bool = True, 
               write_meta: bool = True
               ):
        
        from ..io import dotWrite
        
        with dotWrite(path, append=append) as w:
            w.write(self, write_struct=write_struct, write_meta=write_meta)

    
    def to_fasta(self, path: Union[str, Path], *, 
                 append: bool = False, 
                 ):
        from ..io import fastaWrite
        
        with fastaWrite(path, append=append) as w:
            w.write(self)


    def to_bpseq(self, path: Union[str, Path], *, 
                 write_meta: bool = True
                 ):
        from ..io import bpseqWrite
        
        with bpseqWrite(path) as w:
            w.write(self, write_meta=write_meta)