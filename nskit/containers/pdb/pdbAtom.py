import numpy as np
from ...exceptions import InvalidPDB


    
class PdbAtom:
    __slots__ = ("is_hetatm", "atomn", "name", "altloc", 
                 "mol_name", "chain", "moln", 
                 "insert_code", 
                 "occupancy", "temp", 
                 "segment", 
                 "element", "charge", 
                 "coords")
    
    def __init__(
                self,
                is_hetatm: bool, atomn: int, name: str, altloc: str,
                mol_name: str, chain: str, moln: int,
                insert_code: str,
                x: float, y: float, z: float,
                occupancy: float, temp: float,
                segment: str,
                element: str, charge: int
                ):
        
        self.is_hetatm = is_hetatm
        self.atomn = atomn
        self.name = name
        self.altloc = altloc
        self.mol_name = mol_name
        self.chain = chain
        self.moln = moln
        self.insert_code = insert_code
        self.occupancy = occupancy
        self.temp = temp
        self.segment = segment
        self.element = element
        self.charge = charge
        self.coords = np.array([x,y,z], dtype=np.float32)
    
    
    def __str__(self):
        atom_type = "HETATM" if self.is_hetatm else "ATOM"
        
        if len(self.name)==4:
            name = self.name
        elif len(self.name)==1:
            name = f" {self.name}  "
        elif len(self.element)==1: # H C
            name = f" {self.name:<3}"    
        else: # two characters element (Cl, Fe ...)
            name = self.name.ljust(4, ' ')    
            
        occupancy = f"{self.occupancy:>6.2f}" if isinstance(self.occupancy, float) else " "*6
        temp = f"{self.temp:>6.2f}" if isinstance(self.temp, float) else " "*6
        charge = ['  ', '+ ', '- '][self.charge]
        
        return (f"{atom_type:<6}{self.atomn:>5} "
                f"{name}{self.altloc:>1}"
                f"{self.mol_name:>3} {self.chain}{self.moln:>4}"
                f"{self.insert_code:>1}   "
                f"{self.coords[0]:>8.3f}{self.coords[1]:>8.3f}{self.coords[2]:>8.3f}"
                f"{occupancy}{temp}      "
                f"{self.segment:<4}{self.element:>2}{charge}"
               )
    
    @staticmethod
    def _default_element_derive_func(is_hetatm: bool, name: str, mol_name: str, chain: str):
        return name[0]
    
    @classmethod
    def from_pdb_line(cls,
                      line,
                      derive_element: bool = False,
                      element_derive_func = None
                     ):
        """
        https://www.biostat.jhsph.edu/~iruczins/teaching/260.655/links/pdbformat.pdf
        """
        if derive_element and element_derive_func is None:
            element_derive_func = PdbAtom._default_element_derive_func
            
        line = line.ljust(80, ' ')
        
        is_hetatm = line.startswith("HETATM")  # ATOM or HETATM
        atomn = int(line[6:11].strip())        # Atom serial number
        name = line[12:16].strip()             # Atom name
        altloc = line[16]                      # Alternate location indicator
        mol_name = line[17:20].strip()         # Residue/mol name
        chain = line[21]                       # Chain identifier
        moln = int(line[22:26].strip())        # Residue sequence number
        insert_code = line[26]                 # Code for insertions of residues
        x = float(line[30:38].strip())         # X
        y = float(line[38:46].strip())         # Y
        z = float(line[46:54].strip())         # Z
        
        if (occupancy:=line[54:60].strip()):   # Occupancy
            occupancy = float(occupancy)
            
        if (temp:=line[60:66].strip()):        # Temperature factor
            temp = float(temp)
        
        segment = line[72:76]                  # Segment identifier
        element = line[76:78].strip()          # Element symbol
        if element=="":
            if not derive_element:
                raise InvalidPDB(f"Atom {atomn} {name} has no element field.")
            element = element_derive_func(is_hetatm, name, mol_name, chain)     
        
        charge = line[78:80].strip()           # Charge
        if charge=='': charge = 0
        elif charge=='+': charge = 1
        elif charge=='-': charge = -1
        else:
            raise InvalidPDB(f"Invalid atom charge '{charge}' in atom {atomn} {name}, expected + - or nothing.")
        
        return PdbAtom(is_hetatm, atomn, name, altloc, mol_name, chain, moln, insert_code,
                    x, y, z, occupancy, temp, segment, element, charge)