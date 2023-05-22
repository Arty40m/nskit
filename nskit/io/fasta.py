from typing import Optional, Union, List
from pathlib import Path

from ._DotFasta import DotFastaRead, DotFastaWrite
from ..parse_na import NA
from ..nucleic_acid import NucleicAcid
from ..exceptions import *



class FastaRead(DotFastaRead):
    
    def __init__(self, file: Union[str, Path], *,  
                 raise_na_errors: bool = False,
                 upper_sequence: bool = True,
                ):
        
        super().__init__(file)

        self.raise_na_errors = raise_na_errors
        self.upper_sequence = upper_sequence
    
    
    def _make_na(self, lines: List[str], last_name_idx: int) -> Optional[NucleicAcid]:
        if len(lines)<2:
            raise InvalidFasta(f"Empty structure (structure at line {last_name_idx})")
            
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
    
    
class FastaWrite(DotFastaWrite):
    
    def __init__(self, file, *, 
                 append: bool = False,
                ):
        
        super().__init__(file, append=append)
        
        
    def write(self, data: Union[NucleicAcid, tuple, list]):
        
        if isinstance(data, NucleicAcid):
            name = data.name or 'Seq'
            seq = data.seq
            
        elif isinstance(data, (tuple, list)):
            if len(data)!=2:
                raise ValueError(f"Tuple of structure must contain name and sequence")
                
            name = data[0]
            seq = data[1]
        
        else:
            raise ValueError(f"Data must be NucleicAcid or tuple/list with data, got {type(data)}")
            
        # split seq to chunks of 80 nb
        chunks = []
        i=0
        for i in range(len(seq)//80):
            chunks.append(seq[i*80 : (i+1)*80])

        if len(seq)%80!=0:
            chunks.append(seq[(i+1)*80:])


        self._file.write(f">{name}\n")
        for chunk in chunks:
            self._file.write(f"{chunk}\n")
                
                
                