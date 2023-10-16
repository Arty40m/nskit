from typing import Iterator, Optional, Union
import os
from pathlib import Path
from io import TextIOWrapper

from ..containers import NucleicAcid
from ..exceptions import InvalidStructure



META_SEPARATOR = ": "


class bpseqRead:
    
    def __init__(self, file: Union[str, Path, TextIOWrapper], *, 
                 raise_na_errors: bool = False, 
                 file_as_name: bool = False, 
                ):
        
        if isinstance(file, (str, Path)):
            self._file = open(file)
        elif isinstance(file, TextIOWrapper):
            self._file = file
        else:
            raise TypeError(f"Invalid file type. Accepted - string, Path, TextIOWrapper")
            
        self.raise_na_errors = raise_na_errors
        self.name = os.path.basename(file).split('.')[0] if file_as_name else None
            
        
    def __enter__(self):
        return self

    
    def close(self):
        self._file.close()


    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        
        
    def read(self) -> Optional[NucleicAcid]:
        meta = {}
        seq = []
        pairs = {}
        prev_idx = 0
        
        # parse
        for i, l in enumerate(self._file):
            line = l.strip()
            if not line: continue
            
            if line[0].isnumeric(): # pair
                idx, nb, compl = line.split()
                idx, compl = int(idx), int(compl)
                seq.append(nb)
                
                if (idx-prev_idx)!=1:
                    raise InvalidStructure(f"Nucleic base indexes must be sequential")

                if idx==compl:
                    raise InvalidStructure(f"Nucleic base {idx-1} is self bounded")
                
                prev_idx = idx
                if not compl:
                    continue
                    
                pairs[idx-1] = compl-1
                    
            else: # meta inf
                toks = line.split(META_SEPARATOR)
                if len(toks)!=2: continue
                k, v = toks
                meta[k] = v
                
        # validate
        for o, e in pairs.items():
            # o : e
            # e : o
            if (e<0 or e>=len(seq)) or (o<0 or o>=len(seq)):
                if self.raise_na_errors:
                    raise InvalidStructure(f"Paired nucleic base index out of range")
                return None
            
            expected_o = pairs.get(e)
            if expected_o is None:
                if self.raise_na_errors:
                    raise InvalidStructure(f"Directed bond between nucleic base {o} and {e}")
                return None
            
            if expected_o!=o:
                if self.raise_na_errors:
                    raise InvalidStructure(f"Nucleic base {e} has two bonds")
                return None

        # make na
        na = NucleicAcid()
        if self.name: na.name = self.name
        if meta: na.meta.update(meta)
            
        for i, nb in enumerate(seq):
            _ = na._add_node(nb)
            
        for o, e in pairs.items():
            na._add_bond(o, e)
            
        return na
    
    
class bpseqDirRead:
    
    def __init__(self, Dir: Union[str, Path], *, 
                 raise_na_errors: bool = False, 
                 file_as_name: bool = False, 
                ):
            
        self._dir = Path(Dir)
        
        self.raise_na_errors = raise_na_errors
        self.file_as_name = file_as_name

        self._iterator = self._iterate()
        
        
    def __enter__(self):
        return self
    
    
    def __exit__(self, exc_type, exc_value, traceback):
        ...


    def __len__(self):
        count = 0
        for file in os.listdir(self._dir):
            if file.endswith(".bpseq"):
                count += 1
        return count
    

    def _iterate(self):
        for file in os.listdir(self._dir):
            if not file.endswith(".bpseq"):
                continue
            
            path = self._dir/file
            with bpseqRead(path, 
                           raise_na_errors=self.raise_na_errors, 
                           file_as_name=self.file_as_name, 
                          ) as f:
                yield f.read()
        
    
    def __iter__(self) ->  Iterator[Optional[NucleicAcid]]:
        return self._iterator
    

    def __next__(self) -> Optional[NucleicAcid]:
        return next(self._iterator)
                
                
class bpseqWrite:
    
    def __init__(self, file: Union[str, Path]):
        if isinstance(file, (str, Path)):
            self._file = open(file, 'w')
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
        
        
    def write(self, na: NucleicAcid, *,
              write_meta: bool = True
             ):
        
        if not isinstance(na, NucleicAcid):
            raise TypeError("Can write only NucleicAcid graph")
            
        if write_meta and na.meta:
            for k, v in na.meta.items():
                self._file.write(f"{k}{META_SEPARATOR}{str(v)}\n")
                
        for i, nb in enumerate(na.seq):
            compl = c+1 if (c:=na.complnb(i)) is not None else 0
            self._file.write(f"{i+1} {nb} {compl}\n")
            
            
            