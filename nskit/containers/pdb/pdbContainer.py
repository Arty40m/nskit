from typing import Union, List
import numpy as np
from .pdbAtom import PdbAtom
from .pdbMolecule import PdbMolecule



class PDB:
    def __init__(self, molecules: List[PdbMolecule]):
        ...