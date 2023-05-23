from functools import cached_property
from typing import List, Optional, Tuple
import numpy as np
import numpy

from .graph import Graph



HELIX_ORDER = {0:'()', 1:'[]', 2:'{}', 3:'<>', 4:'Aa', 5:'Bb', 6:'Cc', 7:'Dd', 8:'Ee', 9:'Ff'}


class NucleicAcidGraph(Graph):
    
    def __init__(self):
        super().__init__()
        
        
    def complnb(self, n: int) -> Optional[int]:
        for m, bond_type in self._bonds[n].items():
            if bond_type:
                return m
        return None
    

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
    
    
    @cached_property
    def helixes(self) -> List[tuple]:
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
        
        
    def _helix_intersect(self, current_helix, prev_helix) -> bool:
        # .(((..[[..)))..]].
        return current_helix[0][-1]<prev_helix[1][0] and current_helix[1][0]>prev_helix[1][-1]


    @cached_property
    def helix_orders(self) -> Tuple[int]:
        helixes = self.helixes
        orders = [0 for _ in range(len(helixes))]

        for i, h in enumerate(helixes):
            intersection_orders = []
            for prev_h in range(i):
                if self._helix_intersect(h, helixes[prev_h]):
                    intersection_orders.append(orders[prev_h])

            intersection_orders = set(intersection_orders)
            if len(intersection_orders)>=10:
                break
            lowest_helix_order = min(set(range(10)) - intersection_orders) # lowest type of brackets
            orders[i] = lowest_helix_order

        return tuple(orders)


    def assemble_dot_structure(self) -> str:
        helixes = self.helixes
        orders = self.helix_orders
        dots = ['.' for _ in range(len(self))]

        for i, h in enumerate(helixes):
            brackets = HELIX_ORDER[orders[i]]
            for c in range(len(h[0])):
                dots[h[0][c]] = brackets[0] # open bracket
                dots[h[1][-(1+c)]] = brackets[1] # close bracket

        return ''.join(dots)
    

    @property
    def knots(self) -> Tuple[int]:
        return tuple([i for i in range(len(self.helixes)) if self.helix_orders[i]])
    

    @property
    def knot_helixes(self) -> List[tuple]:
        return [h for i, h in enumerate(self.helixes) if i in self.knots]

    
    @property
    def knot_pairs(self) -> List[tuple]:
        pairs = []
        for h in self.knot_helixes:
            for i in range(len(h[0])):
                pairs.append((h[0][i], h[1][-(1+i)]))

        return pairs
    

    def is_knot(self) -> bool:
        if self.helix_orders:
            return max(self.helix_orders)>0
        return False


    def get_adjacency(self) -> numpy.array:
        slen = len(self)
        adj = np.zeros((slen, slen), dtype=np.int32)
        
        for o, e in self.pairs:
            adj[o, e] = 1
            adj[e, o] = 1
            
        return adj