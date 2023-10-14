from typing import Iterator, Iterable, Optional, Union, List
from pathlib import Path

from .dotLines import dotLinesRead, dotLinesWrite
from ..parse_na import NA
from ..containers import NucleicAcid
from ..exceptions import *



class fastaRead(dotLinesRead):
    
    def __init__(self, file: Union[str, Path], *,  
                 raise_na_errors: bool = False,
                 upper_sequence: bool = False,
                ):
        
        super().__init__(file)

        self.raise_na_errors = raise_na_errors
        self.upper_sequence = upper_sequence
        
        self._na_iterator = self._na_iterate()
        
        
    def _na_iterate(self):
        for i, lines in enumerate(self._iterator):
            yield self._make_na(lines, i)
        
    
    def __iter__(self) -> Iterator[Optional[NucleicAcid]]:
        return self._na_iterator
                    
                    
    def __next__(self) -> Optional[NucleicAcid]:
        return next(self._na_iterator)
    
    
    def _make_na(self, lines: List[str], last_na_idx: int) -> Optional[NucleicAcid]:
        if len(lines)<2:
            raise InvalidFasta(f"Empty structure (structure at line {last_na_idx})")
            
        name = lines[0].strip(">")
        seq = ''.join([l.strip("*") for l in lines[1:] if not l.startswith(";")])
        
        # make NA
        try:
            na = NA(seq, 
                    name=name, 
                    upper_sequence = self.upper_sequence,
                   )
            
        except Exception as e:
            if isinstance(e, (InvalidSequence, InvalidStructure)):
                if self.raise_na_errors:
                    raise e
                return None
            raise e
            
        return na
    
    
class fastaWrite(dotLinesWrite):
    
    def __init__(self, file, *, 
                 append: bool = False,
                ):
        
        super().__init__(file, append=append)
        
        
    def write(self, na: NucleicAcid):
        
        if not isinstance(na, NucleicAcid):
            raise ValueError(f"Data must be NucleicAcid container, got {type(na)}")
            
        name = f">{na.name}" or '>Seq'
        seq = na.seq
        
        # split seq to chunks of 80 nb
        lines = [name]
        i=0
        for i in range(len(seq)//80):
            lines.append(seq[i*80 : (i+1)*80])

        if len(seq)%80!=0:
            lines.append(seq[(i+1)*80:])

        super().write(lines)
                
                
                
                
                
                