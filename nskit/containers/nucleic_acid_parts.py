class Helix:
    def __init__(self, open_chain, close_chain):
        self.__opc = open_chain
        self.__clc = close_chain
    
    
    @property
    def opc(self):
        return self.__opc
    
    
    @property
    def clc(self):
        return self.__clc
    
    
    def __iter__(self):
        return zip(self.__opc, self.__clc[::-1])
    
    
    def __len__(self):
        return len(self.__opc)
    
    
    def __str__(self):
        return (f"-({', '.join([str(i) for i in self.__opc])})\n"
                f" ({', '.join([str(i) for i in self.__clc[::-1]])})-")
    
    
    def __repr__(self):
        r = (f"5'\\{' '*(len(self)*2 - 1)}/3'\n"
             f"   {' '.join([str(i) for i in self.__opc])}\n"
             f"   {' '.join(['|']*len(self))}\n"
             f"   {' '.join(str(i) for i in self.__clc[::-1])}\n"
             f"3'/{' '*(len(self)*2 - 1)}\\5'\n")
        
        return r
    
    
class Loop:
    def __init__(self, nodes):
        self._nodes = nodes
        self._loop = [nodes[0][0]]
        for i, n in enumerate(nodes, start=1):
            if isinstance(n, int):
                self._loop.append(n)
            else:
                self._loop.append(n[0])
                self._loop.append(n[1])
        self._loop.append(nodes[0][1])
                
        self._branches = [n for n in nodes if isinstance(n, tuple)]
        
        
    @property
    def nodes(self):
        return self._nodes
    
    
    @property
    def loop(self):
        return self._loop
    
    
    @property
    def branches(self):
        return self._branches
    
    
    def __len__(self):
        return len(self._loop)
    
    
    def __iter__(self):
        return self._nodes
    
    
    def __str__(self):
        tokens = [""]*(1 + len(self._nodes))
        tokens[0] = f"-{self._loop[0]})"
        tokens[-1] = f"({self._loop[-1]}-"
        
        for i in range(1, len(self._nodes)):
            n = self._nodes[i]
            if isinstance(n, int):
                tokens[i] = str(n)
            else:
                tokens[i] = f"({n[0]}-{n[1]})"

        return ", ".join(tokens)
    
    
    def __repr__(self):
        return str(self)
    
    
class Hairpin(Loop):
    ...
    
    
class InternalLoop(Loop):
    ...
    
    
class Junction(Loop):
    ...
    
    
    
    
    
    
    
    