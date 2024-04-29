from typing import Union, List
import numpy as np
from .pdbAtom import PdbAtom
from .pdbMolecule import PdbMolecule, PdbResidue


    
class PDBCompounds:
    __slots__ = ("__comps", )
    
    def __init__(self):
        self.__comps = []
        
        
    def __getitem__(self, i: int):
        return self.__comps[i]
    
    def __len__(self):
        return len(self.__comps)
    
    def __iter__(self):
        return iter(self.__comps)
        
        
    def add(self, compound: Union[PdbMolecule, PdbResidue, "PDBChain"]):
        if not isinstance(compound, (PdbMolecule, PdbResidue, PDBChain)):
            raise ValueError(f"Can not add {type(compound)} to PDBCompounds, expected PdbMolecule, PdbResidue, PDBChain.")
            
        self.__comps.append(compound)
        
        
class PDBChain(PDBCompounds):
    def __init__(self):
        super().__init__()
        
    
    def add(self, residue: PdbResidue):
        if not isinstance(residue, PdbResidue):
            raise ValueError(f"Can not add {type(residue)} to PDBChain, expected PdbResidue.")
            
        if len(self):
            res = self[0]
            if residue.chain!=res.chain:
                raise ValueError(f"All residues of a chain must have the same chain name ({res.chain}), got {residue.chain}.")
            
        super().add(residue)
        
        
class PDB(PDBCompounds):
    def __init__(self):
        super().__init__()
        
        
class PDBModels:
    __slots__ = ("__models", )
    
    def __init__(self, models: List[PDB]):
        self.__models = models
        
    def __getitem__(self, i: int):
        return self.__models[i]
    
    def __len__(self):
        return len(self.__models)