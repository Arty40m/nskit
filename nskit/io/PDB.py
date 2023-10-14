import numpy as np
import os
from typing import Union, List
from pathlib import Path
from io import TextIOWrapper


from ..containers import NucleicAcid
from ._pdbRead import pdbRead
from ..exceptions import InvalidPDB, InvalidSequence



NA_NAMES = {"A", "U", "G", "C", "T"}

DIRECTION_ATOMS = {
        'A': ('N1', 'N9'),
        'U': ('N3', 'C6'),
        'G': ('N1', 'N9'),
        'C': ('N3', 'C6'),
        'T': ('N3', 'C6'),
    }

DONOR_ACCEPTOR_GROUPS = {
        'A':{
            'donor':(('H61', 'N6'), ('H62', 'N6')),
            'acceptor':('N1', 'N3', 'N7')
        },
        'U':{
            'donor':(('H3', 'N3'),),
            'acceptor':('O2', 'O4')
        },
        'G':{
            'donor':(('H1', 'N1'), ('H21', 'N2'), ('H22', 'N2')),
            'acceptor':('O6', 'N3', 'N7')
        },
        'C':{
            'donor':(('H41', 'N4'), ('H42', 'N4')),
            'acceptor':('O2', 'N3')
        },
        'T':{
            'donor':(('H3', 'N3'),),
            'acceptor':('O2', 'O4')
        },
    }

CLOSE_RESIDUES_THRESHOLD = 25
MIN_DIRECTION_ANGLE = 90
    
MIN_BOND_ANGLE = 120
MAX_BOND_ANGLE = 180

SIGMA10 = 1.9**10
SIGMA12 = 1.9**12
MIN_ENERGY_THRESHOLD = -1.0493852218717947


def norm(x):
    return x/np.linalg.norm(x)


class pdbParse:

    def __init__(self, 
                 close_residues_threshold: float = CLOSE_RESIDUES_THRESHOLD, 
                 min_direction_angle: float = MIN_DIRECTION_ANGLE, 
                 min_bond_angle: float = MIN_BOND_ANGLE, 
                 max_bond_angle: float = MAX_BOND_ANGLE, 
                 min_energy_threshold: float = MIN_ENERGY_THRESHOLD
                 ):
        
        self.close_residues_threshold = close_residues_threshold
        self.min_direction_angle = min_direction_angle
        self.min_bond_angle = min_bond_angle
        self.max_bond_angle = max_bond_angle
        self.min_energy_threshold = min_energy_threshold


    def parse(
            self, file: Union[str, Path, TextIOWrapper], *, 
            assert_non_sequential: bool = False,
            split_chain_by_name: bool = False,
            
            ignore_nonstandard_residues: bool = False, 
            with_energy_matrix: bool = False, 
            return_single_chain: bool = False
              ) -> List[NucleicAcid]:
        
        with pdbRead(file, 
                     assert_non_sequential=assert_non_sequential, 
                     split_chain_by_name=split_chain_by_name) as f:
            chains = f.read()['nas']

        if chains is None:
            raise InvalidPDB(f"PDB file does not have any nucleic acid chain")
        
        nas = []
        for chain in chains:
            name = os.path.basename(file).split('.')[0]

            # sequence
            seq = []
            for res in chain:
                res_name = res.res_name.strip(" D35")
                
                if len(res_name)!=1:
                    if ignore_nonstandard_residues:
                        res_name = res_name[-1]
                    else:
                        raise InvalidSequence(f"Residue with invalid name: {res_name}")
                
                seq.append(res_name)
            seq = ''.join(seq)

            # structure
            energy_matrix = self.get_bond_energies(chain)
            to_adj = energy_matrix.copy() if with_energy_matrix else energy_matrix
            adj = self.min_energy_adjacency(to_adj)
            na = NucleicAcid.from_adjacency(adj, seq=seq, name=name, trust_adj=True)

            if with_energy_matrix:
                nas.append((na, energy_matrix))
            else:
                nas.append(na)

            if return_single_chain:
                break

        return nas
    
    
    def min_energy_adjacency(self, M):
        N = M.shape[0]
        adj = np.zeros((N, N), dtype=np.float32)
        
        while True:
            m = np.argmin(M)
            r = m//N
            c = m%N
            
            if M[r,c]>=self.min_energy_threshold:
                break

            adj[r,c] = 1
            adj[c,r] = 1
            M[r] = 1
            M[c] = 1
            M[:, r] = 1
            M[:, c] = 1
            
        return adj
        
    
    def get_bond_energies(self, chain):
        N = len(chain)
        M = np.ones((N, N), dtype=np.float32)
        
        for i in range(N-2):
            ires = chain[i]
            if ires.res_name.strip(" D35") not in NA_NAMES:
                continue
                
            for j in range(i+2, N):
                jres = chain[j]
                if jres.res_name.strip(" D35") not in NA_NAMES:
                    continue
                
                if not self.close_enough(ires, jres):
                    continue
                    
                if not self.correct_pair_direction(ires, jres):
                    continue
                    
                bonds = self.get_complementary_bonds(ires, jres)
                if len(bonds)==0:
                    continue
                    
                e = self.calculate_bond_energy(bonds)
                M[i, j] = e
                M[j, i] = e
                
        return M
                
    
    def calculate_bond_energy(self, bonds):
        E = 0
        for acceptor, dist in bonds:
            r10 = dist**10
            r12 = r10*(dist**2)
            
            E += 25*SIGMA12/r12 - 30*SIGMA10/r10
        
        return E
        
    
    def get_complementary_bonds(self, r1, r2):
        r1_type = r1.res_name.strip(" D35")[-1]
        r2_type = r2.res_name.strip(" D35")[-1]
        
        r1_acceptors = DONOR_ACCEPTOR_GROUPS[r1_type]['acceptor']
        r2_acceptors = DONOR_ACCEPTOR_GROUPS[r2_type]['acceptor']
        r1_donor_groups = DONOR_ACCEPTOR_GROUPS[r1_type]['donor']
        r2_donor_groups = DONOR_ACCEPTOR_GROUPS[r2_type]['donor']
        
        bonds = [] # (acceptor atom type, distance), ...
        self._add_complementary_bonds(bonds, 
                                     donor_res=r1, 
                                     acceptor_res=r2, 
                                     donor_groups=r1_donor_groups, 
                                     acceptors=r2_acceptors)
        
        self._add_complementary_bonds(bonds, 
                                     donor_res=r2, 
                                     acceptor_res=r1, 
                                     donor_groups=r2_donor_groups, 
                                     acceptors=r1_acceptors)
        
        return bonds
    
    
    def _add_complementary_bonds(self, bonds, donor_res, acceptor_res, donor_groups, acceptors):
        for h, hd in donor_groups: # H atom, H donor atom
            h_vec = donor_res.get_atom_vec(h)
            hd_vec = donor_res.get_atom_vec(hd)
            
            for a in acceptors:
                a_vec = acceptor_res.get_atom_vec(a)
                
                H_donor_vec = norm(hd_vec-h_vec)
                H_acceptor_vec = norm(a_vec-h_vec)
                angle = np.arccos(np.dot(H_donor_vec, H_acceptor_vec))*180/np.pi
                
                if angle>self.max_bond_angle or angle<self.min_bond_angle:
                    continue
                
                dist = np.linalg.norm(a_vec - h_vec)
                bonds.append((a, dist))
        
    
    def close_enough(self, ires, jres):
        c1_i_vec = ires.get_atom_vec("C1'")
        c1_j_vec = jres.get_atom_vec("C1'")
        
        vec = c1_j_vec - c1_i_vec
        d = np.linalg.norm(vec)
        
        if d<=self.close_residues_threshold:
            return True
        
        return False
    
    
    def correct_pair_direction(self, ires, jres):
        ir_type = ires.res_name.strip(" D35")[-1]
        jr_type = jres.res_name.strip(" D35")[-1]
        
        ir_dirdot1, ir_dirdot2 = DIRECTION_ATOMS[ir_type]
        jr_dirdot1, jr_dirdot2 = DIRECTION_ATOMS[jr_type]
        
        ir_vec = ires.get_atom_vec(ir_dirdot1) - ires.get_atom_vec(ir_dirdot2)
        jr_vec = jres.get_atom_vec(jr_dirdot1) - jres.get_atom_vec(jr_dirdot2)
        
        angle = np.arccos(np.dot(norm(ir_vec), norm(jr_vec)))*180/np.pi
        return angle>=self.min_direction_angle