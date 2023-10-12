from typing import Union



class SimplifiedLinearGraph:
    __slots__ = ('_nodes', '_bonds')
    
    def __init__(self):
        self._nodes = []
        self._bonds = {}
        
        
    def _add_node(self, node: Union[str, int, tuple]) -> int:
        self._nodes.append(node)
        return len(self._nodes)-1
    
    
    def _add_bond(self, n: int, m: int):
        self._bonds[n] = m
        self._bonds[m] = n
            
    
    def _remove_bond(self, n: int, m: int):
        _ = self._bonds.pop(m)
        _ = self._bonds.pop(n)
        

class Graph:
    __slots__ = ('_nodes', '_bonds')
    
    def __init__(self):
        self._nodes = {}
        self._bonds = {}
        
        
    def _add_node(self, node: Union[str, int, tuple]) -> int:
        idx = len(self._nodes)
        self._nodes[idx] = node
        self._bonds[idx] = {}
        
        return idx
    
    
    def _add_bond(self, n: int, m: int, bond: int):
        self._bonds[n][m] = bond
        self._bonds[m][n] = bond
        
        
    def _remove_node(self, n: int):
        del self._nodes[n]
        del self._bonds[n]
        
        for node in self._bonds:
            neibs = self._bonds[node]
            _ = neibs.pop(n, None)
            
    
    def _remove_bond(self, n: int, m: int):
        _ = self._bonds[n].pop(m)
        _ = self._bonds[m].pop(n)