from typing import Union, List, Dict
from pathlib import Path
from io import TextIOWrapper
import numpy as np

from ..exceptions import InvalidPDB



NA_NAMES = {"A", "U", "G", "C", "T"}

AMIN_NAMES = {
    'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'ILE', 'LYS', 'LEU', 'MET', 'PRO', 'GLN',
    'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR', 
    'ASH',  'ASN', 
    'HID', 'HIE', 'HIP', 'HIS'
}


class Residue():
    def __init__(self, tokens):
        self.res_name = tokens[0][2]
        self.resn = tokens[0][3]
        self.atoms = tuple([t[1] for t in tokens])
        self.coords = np.array([t[4:] for t in tokens], dtype=np.float32)
        

    def get_idx(self, name):
        for i in range(len(self.atoms)):
            if self.atoms[i] == name:
                return i
        
        raise KeyError(f"No such atom {name} in residue {self.res_name} at number {self.resn}")
    
    
    def get_atom_vec(self, name):
        return self.coords[self.get_idx(name)]
    

class pdbRead:
    def __init__(self, file: Union[str, Path, TextIOWrapper], *, 
                 assert_non_sequential: bool = False, 
                 split_chain_by_name: bool = False
                 ):
        if isinstance(file, (str, Path)):
            self._file = open(file)
        elif isinstance(file, TextIOWrapper):
            self._file = file
        else:
            raise TypeError(f"Invalid file type. Accepted - string, Path, TextIOWrapper")
        
        self.assert_non_sequential = assert_non_sequential
        self.split_chain_by_name = split_chain_by_name
        
        
    def __enter__(self):
        return self

    
    def close(self):
        self._file.close()


    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        
        
    def tokenize(self, line): 
        tokens = (
            int(line[6:11].strip()), # Atom serial number
            line[12:16].strip(), # Atom name
            line[17:20].strip(), # Residue name
            int(line[22:26].strip()), # Residue sequence number
            float(line[30:38].strip()), # X
            float(line[38:46].strip()), # Y
            float(line[46:54].strip()), # Z
        ) 
        
        return tokens
    

    def read(self) -> Dict:
        chains = self.parse_chains()
        chains = self.classify_chains(chains)
        return chains
    

    def get_residue_type(self, res_name: str) -> str:
        if res_name in AMIN_NAMES:
            return 'amins'
        
        if res_name.strip(" D35") in NA_NAMES:
            return 'nas'
        
        return 'ligand'


    def classify_chains(self, chains: List[List[Residue]]) -> Dict:
        classes = {"nas":None, "amins":None, "ligands":None}

        for chain in chains:
            typ = self.get_residue_type(chain[0].res_name)

            if classes[typ] is None: 
                classes[typ] = []
            
            classes[typ].append(chain)

        return classes


    def parse_chains(self) -> List:
        chains = []
        last_res_n = 0
        res_tokens = []
        chain_ress = []
        last_res_type = None

        for line in self._file:

            # Add chain
            if (line.startswith("TER") or \
                line.startswith("MODEL") or \
                line.startswith("ENDMDL")):

                # remaining residue
                if len(res_tokens)>0: 
                    chain_ress.append(Residue(res_tokens))
                    res_tokens = []
                    
                if len(chain_ress)>0:
                    chains.append(chain_ress)
                    chain_ress = []
                    last_res_n = 0

            elif line.startswith('ATOM') or (line.startswith('HETATM')):
                tokens = self.tokenize(line)
                current_res_n = tokens[3]
                current_res_type = self.get_residue_type(tokens[2])

                # implicit chain split by different names
                if current_res_type!=last_res_type and \
                        len(res_tokens)>0 and \
                        self.split_chain_by_name:

                    chain_ress.append(Residue(res_tokens))
                    chains.append(chain_ress)
                    res_tokens = []
                    chain_ress = []
                    last_res_n = current_res_n
                    last_res_type = current_res_type
                    res_tokens.append(tokens)
                    continue

                # new residue in the same chain
                if current_res_n!=last_res_n:
                    if (current_res_n-last_res_n)!=1 and self.assert_non_sequential:
                        raise InvalidPDB((f"Residue numbers must be sequential, got "
                                          "{current_res_n} after {last_res_n}"))

                    if len(res_tokens)>0: 
                        chain_ress.append(Residue(res_tokens))
                        res_tokens = []
                    last_res_n = current_res_n
                    last_res_type = current_res_type

                res_tokens.append(tokens)
            
        if len(res_tokens)>0: 
            chain_ress.append(Residue(res_tokens))
            res_tokens = []

        if len(chain_ress)>0:
            chains.append(chain_ress)

        return chains
