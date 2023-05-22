from functools import cached_property
from typing import List, Optional, Tuple
import numpy as np
import numpy

from .graph import Graph



HELIX_ORDER = {0:'()', 1:'[]', 2:'{}', 3:'<>', 4:'Aa', 5:'Bb', 6:'Cc', 7:'Dd', 8:'Ee', 9:'Ff'}


class NucleicAcidGraph(Graph):
    
    def __init__(self):
        super().__init__()
        
        
    @cached_property
    def pairs(self) -> List[Tuple[int, int]]:
        pairs = []
        visited = set()
        
        for n in self._bonds:
            if n in visited:
                continue
                
            for m, bond_type in self._bonds[n].items():
                if bond_type:
                    pairs.append((n, m))
                    visited.add(m)
                    
        return pairs
    
    
    def complnb(self, n: int) -> Optional[int]:
        for m, bond_type in self._bonds[n].items():
            if bond_type:
                return m
        return None
    
    
    @cached_property
    def helixes(self) -> List[Tuple]:
        helixes = []
        op = []
        close = []

        for o, e in self.pairs:
            if len(op)==0:
                op.append(o)
                close.append(e)
                continue

            if o-op[-1]==1 and close[-1]-e==1:
                op.append(o)
                close.append(e)
            else:
                helixes.append((tuple(op), tuple(close[::-1])))
                op = [o]
                close = [e]

        if len(op)!=0:
            helixes.append((tuple(op), tuple(close[::-1])))
        
        return helixes
        
        
    def assemble_dot_structure(self) -> str:
        helixes = self.helixes
        intersection = [0 for _ in range(len(helixes))]
        dots = ['.' for _ in range(len(self))]

        intersect = lambda h, prev_h: h[0][-1]<prev_h[1][0] and h[1][0]>prev_h[1][-1] # .(((..[[..)))..]].

        for i, h in enumerate(helixes):
            helixes_orders = []
            for prev_h in range(i):
                if intersect(h, helixes[prev_h]):
                    helixes_orders.append(intersection[prev_h])

            helixes_orders = set(helixes_orders)
            if len(helixes_orders)>=10:
                break
            lowest_helix_order = min(set(range(10)) - helixes_orders)# lowest type of brackets
            intersection[i] = lowest_helix_order

            brackets = HELIX_ORDER[lowest_helix_order]
            for c in range(len(h[0])):
                dots[h[0][c]] = brackets[0] # open bracket
                dots[h[1][-(1+c)]] = brackets[1] # close bracket

        return ''.join(dots)
    
    
    def get_adjacency(self) -> numpy.array:
        slen = len(self)
        adj = np.zeros((slen, slen), dtype=np.int32)
        
        for o, e in self.pairs:
            adj[o, e] = 1
            adj[e, o] = 1
            
        return adj