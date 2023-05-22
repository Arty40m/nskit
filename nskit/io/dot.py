from typing import Optional, Union, List
from pathlib import Path

from ._DotFasta import DotFastaRead, DotFastaWrite
from ..parse_na import NA
from ..nucleic_acid import NucleicAcid
from ..exceptions import InvalidFasta, InvalidDotBracket, InvalidSequence, InvalidStructure



META_SEPARATOR = ": "


class DotRead(DotFastaRead):
    
    def __init__(self, file: Union[str, Path], *, 
                 raise_na_errors: bool = False,
                 filter_linear_structures: bool = False, 
                 ignore_unclosed_bonds: bool = False, 
                 upper_sequence: bool = True,
                ):
        
        super().__init__(file)
        
        self.raise_na_errors = raise_na_errors
        self.filter_linear_structures = filter_linear_structures
        self.ignore_unclosed_bonds = ignore_unclosed_bonds
        self.upper_sequence = upper_sequence
    
    
    def _make_na(self, lines: List[str], last_name_idx: int) -> Optional[NucleicAcid]:
        if len(lines)<2:
            raise InvalidFasta(f"Empty structure (structure at line {last_name_idx})")
            
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
                    raise InvalidDotBracket(f"Meta information must contain one separation token '{META_SEPARATOR}' (structure at line {last_name_idx})")
                
                k, v = toks
                meta[k] = v
        
        # make NA
        try:
            na = NA(seq, struct, 
                    name=name, 
                    meta=meta,
                    filter_linear_structures = self.filter_linear_structures, 
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
    
    
class DotWrite(DotFastaWrite):
    
    def __init__(self, file, *, 
                 append: bool = False,
                ):
        
        super().__init__(file, append=append)
        
        
    def write(self, data: Union[NucleicAcid, tuple, list], *, 
              write_struct: bool = True, 
              write_meta: bool = True
             ):
        
        if isinstance(data, NucleicAcid):
            name = data.name or 'Seq'
            seq = data.seq
            struct = data.struct
            meta = data.meta
            
        elif isinstance(data, (tuple, list)):
            if len(data)<2:
                raise ValueError(f"Tuple of structure must contain at least name and sequence")
                
            name = data[0]
            seq = data[1]
            struct = data[2] if len(data)>2 else None
            meta = data[3] if len(data)>3 else None
        
        else:
            raise ValueError(f"Data must be NucleicAcid or tuple/list with data, got {type(data)}")
            
        
        self._file.write(f">{name}\n")
        self._file.write(f"{seq}\n")
        if struct and write_struct:
            self._file.write(f"{struct}\n")
        if meta and write_meta:
            for k, v in meta.items():
                self._file.write(f"{k}{META_SEPARATOR}{str(v)}\n")
                
                
                