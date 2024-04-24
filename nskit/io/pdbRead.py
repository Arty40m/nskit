from typing import Union, List, Dict
from pathlib import Path
from io import TextIOWrapper
import numpy as np

from ..containers.pdb.pdbAtom import PdbAtom
from ..containers.pdb.pdbMolecule import PdbMolecule
from ..containers.pdb.pdbContainer import PDB
from ..exceptions import InvalidPDB



class pdbRead:
    def __init__(self, file: Union[str, Path, TextIOWrapper]):
        if isinstance(file, (str, Path)):
            self._file = open(file)
        elif isinstance(file, TextIOWrapper):
            self._file = file
        else:
            raise TypeError(f"Invalid file type. Accepted - string, Path, TextIOWrapper")
        
    def __enter__(self):
        return self
    
    def close(self):
        self._file.close()

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        
        
    def read(self):
        lines = self._file.readlines()
        lines = list(map(lambda l: l.rstrip(), lines))
        
        # Header
        for i, l in enumerate(lines):
            if l.startswith("ATOM") or \
                l.startswith("HETATM") or \
                 l.startswith("MODEL"):
                break
                
        header = "" if i==0 else "\n".join(lines[:i])
        lines = lines[i:]
            
        # Parse atoms
        mols = []
        models = [0]
        chains = [0]
        mol_atoms = []
        cur_mol_idx = -1
        reset_mol_idx = True
        
        tokens = self.parse_lines(lines)
        for i, t in enumerate(tokens):
            if isinstance(t, PdbAtom):
                atom = t
                
                if reset_mol_idx:
                    cur_mol_idx = atom.moln
                    reset_mol_idx = False
                    
                if atom.moln!=cur_mol_idx:
                    mols.append(self.make_mol(mol_atoms))
                    mol_atoms = [atom]
                    cur_mol_idx = atom.moln
                    continue
                    
                mol_atoms.append(atom)
            
            elif t=="MODEL" and len(mols)>0:
                mols.append(self.make_mol(mol_atoms))
                mol_atoms = []
                models.append(len(mols))
                reset_mol_idx = True
                continue
                
            elif t=="TER" and len(mols)>0:
                mols.append(self.make_mol(mol_atoms))
                mol_atoms = []
                chains.append(len(mols))
                reset_mol_idx = True
                continue
                
        return mols
        
        
    def parse_lines(self, lines):
        tokens = []
        for l in lines:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                tokens.append(PdbAtom.from_pdb_line(l))
            else:
                tokens.append(l)
        
        return tokens
    
    def make_mol(self, atoms):
        m = PdbMolecule()
        for atom in atoms:
            m.add_atom(atom)
        return m