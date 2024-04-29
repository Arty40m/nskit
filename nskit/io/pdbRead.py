from typing import Union, List, Dict
from pathlib import Path
from io import TextIOWrapper
import numpy as np

from ..containers.pdb.pdbAtom import PdbAtom
from ..containers.pdb.pdbMolecule import PdbMolecule, PdbResidue
from ..containers.pdb.pdbContainer import PDB, PDBChain, PDBModels
from ..exceptions import InvalidPDB



NA_NAMES = {"A", "U", "G", "C", "I", 
            "DA", "DT", "DG", "DC"}

AMINOACID_NAMES = {'ALA', 'CYS', 'ASP', 'GLU', 
                   'PHE', 'GLY', 'ILE', 'LYS', 
                   'LEU', 'MET', 'PRO', 'GLN', 
                   'ARG', 'SER', 'THR', 'VAL', 
                   'TRP', 'TYR', 'ASH',  'ASN', 
                   'HID', 'HIE', 'HIP', 'HIS'}


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
        tokens = self.parse_atoms(lines)
        models = self.split_models(tokens)
        
        for i, model in enumerate(models):
            model = self.parse_mols(model)
            model = self.parse_chains(model)
            pdb = PDB()
            for component in model: pdb.add(component)
            models[i] = pdb
            
        return PDBModels(models)
    
        
    def parse_atoms(self, lines):
        tokens = []
        
        for l in lines:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                atom = PdbAtom.from_pdb_line(l)
                if atom.altloc not in (' ', 'A'): # altloc != ' ' or 'A'
                    continue
                tokens.append(atom)
                
            elif l.startswith("TER") or l.startswith("MODEL") or l.startswith("ENDMDL"):
                tokens.append(l.strip())
                
        return tokens
    
    
    def split_models(self, atoms):
        models = []
        model_state = 0 # 0 - undefined, 1 - opened (after MODEL), 2 - closed (ENDMDL)
        model_components = []
        
        for a in atoms:
            if isinstance(a, PdbAtom) or a.startswith("TER"):
                model_components.append(a)
                
            elif a.startswith("MODEL"):
                if model_state==1:
                    raise InvalidPDB(f"PDB contains second MODEL line ({a}) without closing ENDMDL before.")
                    
                if len(model_components):
                    raise InvalidPDB(f"PDB contains atom or TER lines before opening ({a}) line.")
                    
                model_state = 1
                
            else: # ENDMDL
                if model_state!=1:
                    raise InvalidPDB(f"PDB contains closing ({a}) line without opening MODEL line before.")
                    
                if len(model_components)==0:
                    raise InvalidPDB(f"PDB contains empty model. Closed at line ({a}).")
                    
                models.append(model_components)
                model_components = []
                model_state = 2
                
        if len(model_components):
            if model_state!=0:
                raise InvalidPDB(f"The last PDB model is not closed by ENDMDL.")
            else: # signle model without MODEL and ENDMDL in pdb file
                models.append(model_components)
        
        return models
    
    
    def parse_mols(self, atom_tokens):
        mol_tokens = []
        mol_atoms = []
        reset_mol_idx = True
        
        for at in atom_tokens:
            if isinstance(at, PdbAtom):
                if reset_mol_idx:
                    cur_mol_idx = at.moln
                    reset_mol_idx = False
                    
                if at.moln!=cur_mol_idx:
                    if len(mol_atoms)>0:
                        mol_tokens.append(self.make_mol(mol_atoms))
                    mol_atoms = []
                    cur_mol_idx = at.moln
                    
                mol_atoms.append(at)
                
            else:
                if len(mol_atoms)>0:
                    mol_tokens.append(self.make_mol(mol_atoms))
                mol_atoms = []
                reset_mol_idx = True
                mol_tokens.append(at)
                
        return mol_tokens
                
        
    def parse_chains(self, mol_tokens):
        chains = []
        chain_mols = []
        reset_chain = True
        
        for m in mol_tokens:
            if isinstance(m, PdbResidue):
                if reset_chain:
                    cur_chain = m.chain
                    reset_chain = False
                    
                if m.chain!=cur_chain:
                    chains.append(self.make_chain(chain_mols))
                    chain_mols = []
                    cur_chain = at.chain
                    
                chain_mols.append(m)
                
            elif isinstance(m, PdbMolecule):
                if len(chain_mols):
                    chains.append(self.make_chain(chain_mols))
                    chain_mols = []
                    reset_chain = True
                chains.append(m)
                
            else: # TER
                if len(chain_mols)==0:
                    raise InvalidPDB(f"PDB contains empty chain at ({m}).")
                    
                chains.append(self.make_chain(chain_mols))
                chain_mols = []
                reset_chain = True
            
        return chains
        
        
    def make_mol(self, atoms):
        if (atoms[0].mol_name in NA_NAMES) or (atoms[0].mol_name in AMINOACID_NAMES):
            m = PdbResidue()
        else:
            m = PdbMolecule()
            
        for atom in atoms:
            m.add_atom(atom)
        return m
    
    
    def make_chain(self, mols):
        chain = PDBChain()
        for m in mols:
            chain.add(m)
        
        return chain