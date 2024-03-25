class Helix:
    
    __slots__ = '__opc', '__clc'
    
    def __init__(self, open_chain, close_chain):
        self.__opc = tuple(open_chain)
        self.__clc = tuple(close_chain)
    
    
    @property
    def opc(self):
        return self.__opc
    
    
    @property
    def clc(self):
        return self.__clc
    
    
    @property
    def root(self):
        return self[0]
    
    
    def __iter__(self):
        return zip(self.__opc, self.__clc[::-1])
    
    
    def __len__(self):
        return len(self.__opc)
    
    
    def __getitem__(self, key): 
        return self.opc[key], self.clc[-(key+1)]
    
    
    def __str__(self):
        return "He-"+str(tuple(iter(self)))
    
    
    def __repr__(self):
        return str(self)
    

def _make_loop(nodes, knots=None):
    c = len([1 for n in nodes if isinstance(n, tuple)])
    if c==1:
        return Hairpin(nodes, knots)
    
    elif c>2:
        return Junction(nodes, knots)
    
    elif isinstance(nodes[1], tuple) or isinstance(nodes[-1], tuple):
        return Bulge(nodes, knots)
    
    else:
        return InternalLoop(nodes, knots)
        
    
class Loop:
    
    __slots__ = '__nodes', '__knots'
    
    def __init__(self, nodes, knots=None):
        self.__nodes = nodes
        self.__knots = knots
        
        
    @property
    def nodes(self):
        return self.__nodes
    
    
    @property
    def knots(self):
        return self.__knots
    
    
    def has_knot(self):
        return self.__knots is not None
    
    
    @property
    def nts(self):
        nts = [self.__nodes[0][0]]
        for i in range(1, len(self.__nodes)):
            n = self.__nodes[i]
            if isinstance(n, int):
                nts.append(n)
            else:
                nts.append(n[0])
                nts.append(n[1])
        nts.append(self.__nodes[0][1])
        return nts
    
    
    @property
    def branches(self):
        return tuple((n for n in self.__nodes if isinstance(n, tuple)))
    
    
    @property
    def root(self):
        return self.__nodes[0]
    
    
    def __len__(self):
        return len(self.__nodes) + len(self.branches)
    
    
    def __iter__(self):
        return self.nts
    
    
    def __str__(self):
        return str(self.__nodes)
    
    
    def __repr__(self):
        return str(self)
    
    
class Hairpin(Loop):
    
    __slots__ = tuple()
    
    def __init__(self, nodes, knots):
        super().__init__(nodes, knots)
        
        
    def __str__(self):
        return "Hp-"+super().__str__()
    
    
class InternalLoop(Loop):
    
    __slots__ = "short_side", "long_side"
    
    def __init__(self, nodes, knots):
        super().__init__(nodes, knots)
        
        a = 1
        for i in range(2, len(self.nodes)):
            if isinstance(self.nodes[i], tuple):
                break
            a+=1
        
        b = len(self.nodes) - 2 - a
        self.short_side = min(a, b)
        self.long_side = max(a, b)
        
        
    def is_symmetric(self):
        return self.short_side==self.long_side
    
    
    def __str__(self):
        return "I-"+super().__str__()
        
        
class Bulge(Loop):
    
    __slots__ = tuple()
    
    def __init__(self, nodes, knots):
        super().__init__(nodes, knots)
        
        
    def __str__(self):
        return "B-"+super().__str__()
        
        
class Junction(Loop):
    
    __slots__ = tuple()
    
    def __init__(self, nodes, knots):
        super().__init__(nodes, knots)
        
        
    def __str__(self):
        return "J-"+super().__str__()
        
        
        
        
        
    
    
    