from typing import Iterator, Iterable, Optional, Union, List
from pathlib import Path
from io import TextIOWrapper

from ..containers import NucleicAcid
from ..exceptions import InvalidFasta



class dotLinesRead:

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
        i = 0
        
        while True:
            line = self._file.readline()
            if not line: break # EOF
            
            i += 1
            line = line.strip()
            if not line: continue # empty line
            
            if line.startswith(">"): # new na
                yield tuple(lines)
                lines = [line]
                continue
            
            lines.append(line)
            
        yield tuple(lines)

        
    def __iter__(self) -> Iterator[Iterable[str]]:
        return self._iterator
                    
                    
    def __next__(self) -> Iterable[str]:
        return next(self._iterator)
            

class dotLinesWrite:
    
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


    def write(self, data: Iterable[str], **kwargs):
        if len(data)<2:
            raise ValueError(f"At least two lines required")
            
        if not all([isinstance(d, str) for d in data]):
            raise ValueError(f"All passed data must be strings")
            
        for i in range(0, len(data)):
            self._file.write(f"{data[i]}\n")
        
        
        
        
        
        
        
        
        
        