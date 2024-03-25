from functools import cached_property
from typing import List, Optional, Tuple, Union
import numpy as np
import numpy

from .graph import SimplifiedLinearGraph
from .nucleic_acid_fragments import Helix, _make_loop, Hairpin, InternalLoop, Bulge, Junction



HELIX_ORDER = {0:'()', 1:'[]', 2:'{}', 3:'<>', 4:'Aa', 5:'Bb', 6:'Cc', 7:'Dd', 8:'Ee', 9:'Ff'}


class NucleicAcidGraph(SimplifiedLinearGraph):

    GRAPH_CACHE_KEYS = ('struct', 'pairs', 'helixes', 'helix_orders', 'knots', 'knot_helixes', 'knot_pairs', 'loops')
    
    def __init__(self):
        super().__init__()
        
        
    def complnb(self, n: int) -> Optional[int]:
        n = n+len(self) if n<0 else n
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
            
        if self.complnb(n):
            raise ValueError(f"Nb {n} already has complementary bond")
        if self.complnb(m):
            raise ValueError(f"Nb {m} already has complementary bond")

        self._add_bond(n, m)
        self.clear_graph_cache()
    

    def split(self, n: int, m: int, clear_cache: bool = True):
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
        if clear_cache:
            self.clear_graph_cache()
        
        
    def fix_sharp_hairpins(self, min_pin_size: int = 1):
        """
        Breaks complementary bonds at the end of each helix until all hairpins are greater or equal then min_pin_size.
        """
        if not isinstance(min_pin_size, int) or (min_pin_size<1):
            raise ValueError("min_pin_size must be integer >= 1")
            
        check_helixes=True
        while check_helixes:
            helixes = self.helixes
            orders = self.helix_orders
            check_helixes = False
            
            for i, h in enumerate(helixes):
                if orders[i]: # pseudoknot helix
                    continue

                pin_size = h.clc[0] - h.opc[-1] - 1
                if pin_size>=min_pin_size:
                    continue

                n_pairs_to_remove = (dif:=min_pin_size-pin_size)//2 + dif%2
                for j in range(1, n_pairs_to_remove+1):
                    if j>len(h): # end of helix
                        break

                    self.split(*h[-j], clear_cache=False)
                    check_helixes=True
                
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
            
            root = h[-1]
            end_idx = root[1]
            loop = [h[-1],]
            idx = root[0] + 1

            knot_helix = []
            last_knot_idx = None
            last_compl_idx = None
            knot_helixes = []
            while True:
                if idx==end_idx: # end of the loop
                    break

                if ((cidx:=self.complnb(idx)) is not None): # nb has complementary bond
                    # normal helix
                    if idx not in knot_nbs: 
                        loop.append((idx, cidx))
                        idx = cidx+1
                        continue

                    # knot helix
                    if last_knot_idx is None: # first knot nb
                        knot_helix.append(idx)

                    elif abs(last_knot_idx-idx)>1 or abs(last_compl_idx-cidx)>1: # new knot helix
                        knot_helixes.append(tuple(knot_helix))
                        knot_helix = [idx]

                    else: # continue knot
                        knot_helix.append(idx)

                    last_knot_idx = idx
                    last_compl_idx = cidx

                loop.append(idx)
                idx+=1

            if len(knot_helix)>0:
                knot_helixes.append(tuple(knot_helix))
                loops.append(_make_loop(tuple(loop), tuple(knot_helixes)))
            else:
                loops.append(_make_loop(tuple(loop)))

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
    
    
    @cached_property
    def dangling_ends(self) -> Tuple[Tuple[int], Tuple[int]]:
        end5 = []
        end3 = []
        for i in range(len(self)):
            if self.complnb(i) is None:
                end5.append(i)
            else:
                break
                
        for j in range(len(self)-1, i, -1):
            if self.complnb(j) is None:
                end3.append(j)
            else:
                break
                
        return tuple(end5), tuple(end3[::-1])


    def get_adjacency(self) -> numpy.array:
        slen = len(self)
        adj = np.zeros((slen, slen), dtype=np.int32)
        
        for o, e in self.pairs:
            adj[o, e] = 1
            adj[e, o] = 1
            
        return adj