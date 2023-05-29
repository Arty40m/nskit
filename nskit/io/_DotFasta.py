from abc import ABC, abstractmethod
from typing import Iterator, Optional, Union, List
from pathlib import Path
from io import TextIOWrapper

from ..nucleic_acid import NucleicAcid
from ..exceptions import InvalidFasta



class DotFastaRead(ABC):

    def __init__(self, file: Union[str, Path, TextIOWrapper]):
        if isinstance(file, (str, Path)):
            self._file = open(file)
        elif isinstance(file, TextIOWrapper):
            self._file = file
        else:
            raise TypeError(f"Invalid file type. Accepted - string, Path, TextIOWrapper")
        
        self._iterator = self._iterate()
        
        
    def __enter__(self):
        return self

    
    def close(self):
        self._file.close()


    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        
        
    def __len__(self):
        tell = self._file.tell()
        self._file.seek(0)
        count = 0
        for l in self._file:
            if l.startswith(">"):
                count+=1
                
        self._file.seek(tell)
        return count
        
        
    def _iterate(self):
        line = self._file.readline().strip()
        if not line.startswith(">"):
            raise InvalidFasta(f"First line name without '>'")
        
        lines = [line]
        last_name_idx = 1
        i = 0
        
        while True:
            line = self._file.readline()
            if not line: break # EOF
            
            i += 1
            line = line.strip()
            if not line: continue # empty line
            
            if line.startswith(">"): # new na
                yield self._make_na(lines, last_name_idx)
                lines = [line]
                last_name_idx = i+1
                continue
            
            lines.append(line)
            
        yield self._make_na(lines, last_name_idx)

        
    def __iter__(self) -> Iterator[Optional[NucleicAcid]]:
        return self._iterator
                    
                    
    def __next__(self) -> Optional[NucleicAcid]:
        return next(self._iterator)
    
    
    @abstractmethod
    def _make_na(self, lines: List[str], last_name_idx: int) -> Optional[NucleicAcid]:
        """
        Makes NucleicAcid from lines
        """
            

class DotFastaWrite(ABC):
    
    def __init__(self, file: Union[str, Path, TextIOWrapper], *, 
                 append: bool = False,
                ):
        
        if isinstance(file, (str, Path)):
            self._file = open(file, 'a' if append else 'w')
        elif isinstance(file, TextIOWrapper):
            self._file = file
        else:
            raise TypeError(f"Invalid file type. Accepted - string, Path, TextIOWrapper")
        
        
    def __enter__(self):
        return self

    
    def close(self):
        self._file.close()


    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


    @abstractmethod
    def write(self, data: Union[NucleicAcid, tuple, list], **kwargs):
        """
        Writes na to file 
        """