from functools import cached_property
import warnings
from typing import List, Optional, Tuple, Union
import numpy as np
import numpy

from .graph import SimplifiedLinearGraph
from .nucleic_acid_parts import Helix, _make_loop, Hairpin, InternalLoop, Bulge, Junction



HELIX_ORDER = {0:'()', 1:'[]', 2:'{}', 3:'<>', 4:'Aa', 5:'Bb', 6:'Cc', 7:'Dd', 8:'Ee', 9:'Ff'}


class NucleicAcidGraph(SimplifiedLinearGraph):

    GRAPH_CACHE_KEYS = ('struct', 'pairs', 'helixes', 'helix_orders', 'knots', 'knot_helixes', 'knot_pairs', 'loops')
    
    def __init__(self):
        super().__init__()
        
        
    def complnb(self, n: int) -> Optional[int]:
        return self._bonds.get(n, None)
        

    def join(self, n: int, m: int):
        """
        Creates complementary bond between specified nbs.
        """
        if n<0: n = len(self) + n
        if m<0: m = len(self) + m

        if (not 0<=n<len(self)) or (not 0<=m<len(self)):
            raise IndexError(f"Nb index is out of range")

        if n==m: 
            raise ValueError("Can not join nb with itself")
            
        if abs(n-m)==1:
            warnings.warn(f"Nbs {n} and {m} form sharp helix")

        if self.complnb(n):
            raise ValueError(f"Nb {n} already has complementary bond")
        if self.complnb(m):
            raise ValueError(f"Nb {m} already has complementary bond")

        self._add_bond(n, m)
        self.clear_graph_cache()
    

    def split(self, n: int, m: int):
        """
        Breaks complementary bond between specified nbs if exists.
        """
        if n<0: n = len(self) + n
        if m<0: m = len(self) + m
        
        if (not 0<=n<len(self)) or (not 0<=m<len(self)):
            raise IndexError(f"Nb index is out of range")

        if self.complnb(n)!=m:
            raise ValueError(f"Nbs {n} and {m} does not have complementary bond")

        self._remove_bond(n, m)
        self.clear_graph_cache()


    def clear_graph_cache(self):
        for key in self.GRAPH_CACHE_KEYS:
            if key in self.__dict__:
                del self.__dict__[key]
    

    @cached_property
    def pairs(self) -> Tuple[Tuple[int, int]]:
        pairs = []
        used = set()
        
        bond_pairs = sorted(self._bonds.items(), key=lambda x: x[0])
        for o, e in bond_pairs:
            if o in used:
                continue
                
            pairs.append((o, e))
            used.add(e)
        
        return tuple(pairs)
        
    
    @cached_property
    def helixes(self) -> Tuple[Helix]:
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
                helixes.append(Helix(tuple(op), tuple(close[::-1])))
                op = [o]
                close = [e]

        if len(op)!=0:
            helixes.append(Helix(tuple(op), tuple(close[::-1])))
        
        return tuple(helixes)
        
        
    def _helix_intersect(self, current_helix, prev_helix) -> bool:
        # .(((..[[..)))..]].
        return current_helix.opc[-1]<prev_helix.clc[0] and current_helix.clc[0]>prev_helix.clc[-1]


    @cached_property
    def helix_orders(self) -> Tuple[int]:
        helixes = self.helixes
        orders = [0 for _ in range(len(helixes))]
        order_range_set = set(range(len(HELIX_ORDER)))

        for i, h in enumerate(helixes):
            intersection_orders = set()
            for prev_h in range(i):
                if self._helix_intersect(h, helixes[prev_h]):
                    intersection_orders.add(orders[prev_h])

            if len(intersection_orders)>=len(HELIX_ORDER):
                break

            lowest_helix_order = min(order_range_set - intersection_orders) # lowest type of brackets
            orders[i] = lowest_helix_order

        return tuple(orders)


    def assemble_dot_structure(self) -> str:
        helixes = self.helixes
        orders = self.helix_orders
        dots = ['.' for _ in range(len(self))]

        for i, h in enumerate(helixes):
            brackets = HELIX_ORDER[orders[i]]
            for o, e in h:
                dots[o] = brackets[0] # open bracket
                dots[e] = brackets[1] # close bracket

        return ''.join(dots)
    

    @property
    def knots(self) -> Tuple[int]:
        return tuple([i for i in range(len(self.helixes)) if self.helix_orders[i]])
    

    @property
    def knot_helixes(self) -> Tuple[Tuple[int]]:
        return tuple([h for i, h in enumerate(self.helixes) if i in self.knots])

    
    @property
    def knot_pairs(self) -> Tuple[tuple]:
        pairs = []
        for h in self.knot_helixes:
            for p in h:
                pairs.append(p)

        return tuple(pairs)
    

    def is_knot(self) -> bool:
        if self.helix_orders:
            return max(self.helix_orders)>0
        return False
    
    ################################  Loops
    
    @cached_property
    def loops(self) -> Tuple[Union[Hairpin, InternalLoop, Bulge, Junction]]:
        knots = set(self.knots)
        knot_nbs = set()
        for i, j in self.knot_pairs:
            knot_nbs.add(i)
            knot_nbs.add(j)

        loops = []
        for i, h in enumerate(self.helixes):
            if i in knots:
                continue

            end_idx = h.clc[0]
            loop = [h[-1],]
            idx = h.opc[-1] + 1
            while True:
                if idx==end_idx: # end of the loop
                    break

                if ((cidx:=self.complnb(idx)) is not None) and (idx not in knot_nbs):
                    loop.append((idx, cidx))
                    idx = cidx+1
                else:
                    loop.append(idx)
                    idx+=1

            loops.append(_make_loop(loop))
            
        return tuple(loops)
    
    
    @property
    def hairpins(self) -> Tuple[Hairpin]:
        return tuple([l for l in self.loops if isinstance(l, Hairpin)])
    
    
    @property
    def internal_loops(self) -> Tuple[InternalLoop]:
        return tuple([l for l in self.loops if isinstance(l, InternalLoop)])
    
    
    @property
    def bulges(self) -> Tuple[Bulge]:
        return tuple([l for l in self.loops if isinstance(l, Bulge)])
    
    
    @property
    def junctions(self) -> Tuple[Junction]:
        return tuple([l for l in self.loops if isinstance(l, Junction)])


    def get_adjacency(self) -> numpy.array:
        slen = len(self)
        adj = np.zeros((slen, slen), dtype=np.int32)
        
        for o, e in self.pairs:
            adj[o, e] = 1
            adj[e, o] = 1
            
        return adj