from typing import Iterator, Iterable, Optional, Union, List
from pathlib import Path

from .dotLines import dotLinesRead, dotLinesWrite
from ..parse_na import NA
from ..containers.nucleic_acid import NucleicAcid
from ..exceptions import InvalidFasta, InvalidDotBracket, InvalidSequence, InvalidStructure



META_SEPARATOR = ": "


class dotRead(dotLinesRead):
    
    def __init__(self, file: Union[str, Path], *, 
                 raise_na_errors: bool = False, 
                 allow_sharp_helixes: bool = False, 
                 fix_sharp_helixes: bool = False, 
                 ignore_unclosed_bonds: bool = False, 
                 upper_sequence: bool = False,
                ):
        
        super().__init__(file)
        
        self.raise_na_errors = raise_na_errors
        self.allow_sharp_helixes = allow_sharp_helixes
        self.fix_sharp_helixes = fix_sharp_helixes
        self.ignore_unclosed_bonds = ignore_unclosed_bonds
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
            raise InvalidFasta(f"Empty structure at index {last_na_idx}")
            
        name = lines[0].strip(" >")
        seq = lines[1]
        
        if len(lines)==2:
            struct = None
            lines = []
        
        elif META_SEPARATOR in lines[2]:
            struct = None
            lines = lines[2:]
            
        else:
            struct = lines[2].strip()
            lines = lines[3:]
            
        # parse meta
        if not len(lines):
            meta = None
        else:
            meta = {}
            
            for ml in lines:
                toks = ml.strip().split(META_SEPARATOR)
                if len(toks)!=2:
                    raise InvalidDotBracket(f"Meta information must contain one separation token '{META_SEPARATOR}' " 
                                            f"(structure at index {last_na_idx}")
                
                k, v = toks
                meta[k] = v
        
        # make NA
        try:
            na = NA(seq, struct, 
                    name=name, 
                    meta=meta, 
                    allow_sharp_helixes = self.allow_sharp_helixes, 
                    fix_sharp_helixes = self.fix_sharp_helixes, 
                    ignore_unclosed_bonds = self.ignore_unclosed_bonds, 
                    upper_sequence = self.upper_sequence,
                   )
            
        except Exception as e:
            if isinstance(e, (InvalidSequence, InvalidStructure)):
                if self.raise_na_errors:
                    raise e
                return None
            raise e
            
        return na
    
    
class dotWrite(dotLinesWrite):
    
    def __init__(self, file, *, 
                 append: bool = False,
                ):
        
        super().__init__(file, append=append)
        
        
    def write(self, na: NucleicAcid, *, 
              write_struct: bool = True, 
              write_meta: bool = True
             ):
        
        if not isinstance(na, NucleicAcid):
            raise ValueError(f"Data must be NucleicAcid container, got {type(na)}")
        
        name = f">{na.name}" or '>Seq'
        lines = [name, na.seq]
        if na.struct is not None:
            lines.append(na.struct)
        
        if na.meta and write_meta:
            for k, v in na.meta.items():
                lines.append(f"{k}{META_SEPARATOR}{str(v)}\n")
        
        super().write(lines)
                
                
                
                
                