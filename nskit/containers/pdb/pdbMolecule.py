from typing import Union, List
import numpy as np
from .pdbAtom import PdbAtom



class PdbMolecule:
    __slots__ = ("__atoms", "__name_idx_map")
    
    def __init__(self):
        self.__atoms = []
        self.__name_idx_map = {}
        
    def _remap(self):
        self.__name_idx_map.clear()
        self.__name_idx_map = {atom.name:i for i, atom in enumerate(self.__atoms)}
        
        
    def __getitem__(self, i: int):
        return self.__atoms[i]
    
    def __len__(self):
        return len(self.__atoms)
    
    def __iter__(self):
        return iter(self.__atoms)
    
    def __str__(self):
        return "\n".join([str(a) for a in self.__atoms])
    
        
    def add_atom(self, atom: PdbAtom):
        if len(self.__atoms):
            if self.__name_idx_map.get(atom.name) is not None:
                raise ValueError(f"Atom ({atom.atomn}) with name {atom.name} already exists in molecule.")
            
            a = self.__atoms[0]
            if atom.mol_name!=a.mol_name:
                raise ValueError(f"All atoms of a molecule (number {a.moln}) must have the same molecule name ({a.mol_name}), got {atom.mol_name}.")

            if atom.chain!=a.chain:
                raise ValueError(f"All atoms of a molecule (number {a.moln}) must have the same chain name ({a.chain}), got {atom.chain}.")
            
        self.__atoms.append(atom)
        self.__name_idx_map[atom.name] = len(self.__atoms)
        
    def get_atom_idx(self, name: str):
        return self.__name_idx_map.get(name)
    
    def get_atom(self, i: Union[int, str]):
        if isinstance(i, int):
            return self.__atoms[i]
        elif isinstance(i, str):
            return self.__atoms[self.__name_idx_map[i]]
        else:
            raise ValueError(f"Invalid argument of type {type(i)}, accepted: int index or str name.")
            
    def delete_atom(self, i: Union[int, str]):
        if isinstance(i, int):
            return self.__atoms.pop(i)
        elif isinstance(i, str):
            return self.__atoms.pop(self.__name_idx_map[i])
        else:
            raise ValueError(f"Invalid argument of type {type(i)}, accepted: int index or str name.")
        self._remap()
        
    
    @property
    def moln(self):
        return self.__atoms[0].moln
    
    @property
    def chain(self):
        return self.__atoms[0].chain
    
    @property
    def coords(self):
        return np.stack([a.coords for a in self.__atoms])
    
    @coords.setter
    def coords(self, c: np.ndarray):
        if len(self)!=c.shape[0] or c.shape[1]!=3:
            raise ValueError(f"Coords matrix must have shape: ({len(self)}, 3), got {c.shape}")
            
        for i in range(c.shape[0]):
            self.__atoms[i].coords = c[i]
            
            
class PdbResidue(PdbMolecule):
    def __init__(self):
        super().__init__()