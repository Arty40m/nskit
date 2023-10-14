from typing import Union, Iterator
from pathlib import Path
from io import BufferedWriter, BufferedRandom
import math

from ..containers import NucleicAcid
from ..exceptions import InvalidSequence



SUPPORTED_NB_SYMBOLS = set('AGCUTIN')
META_SEPARATOR = "?"
NB_DICT = {'N':0, 'A':1, 'U':2, 'G':3, 'C':4, 'T':5, 'I':6}
INV_NB_DICT = {0:'N', 1:'A', 2:'U', 3:'G', 4:'C', 5:'T', 6:'I'}

format_doc = \
"""
Bytes Nucleic Acid - memory efficient nucleic acid file format.

Maximum supported sequence length - 4095 nbs.
Maximum supported name length - 4095 ascii characters. 
Meta information is added to name in json-like format.
Maximum memory consumption is [1.25*N + 5] bytes for N nbs.

12 bit format specification:

    16 bit - size of current structure in bytes (including these 2 bytes)
    12 bit - na name length (with meta information json)
    12 bit - sequence length

    Name block:
    8 bit - for every ascii character in name

    Sequence block:
    1 bit - paired flag, 1 if nb has a complementary bond
    3 bit - for every nb, encoding nb type
        000 - N    001 - A    010 - U    011 - G
        100 - C    101 - T    110 - I
    Last 4 bit in byte are padding in case of odd number of nbs.

    Pairs block:
    12 bit - index of complementary nb from 3'-end
    3'-end indexes are paired sequentially with 5'-end nbs with complementary falg bit.
    Last 4 bit in byte are padding in case of odd number of complementary pairs.
"""


class bnaWrite:
    
    def __init__(self, file: Union[str, Path, BufferedWriter, BufferedRandom], *, 
                 append: bool = False,
                ):
        
        if isinstance(file, (str, Path)):
            self._file = open(file, 'ab' if append else 'wb')
        elif isinstance(file, (BufferedWriter, BufferedRandom)):
            self._file = file
        else:
            raise TypeError(f"Invalid file type. Accepted - string, Path, BufferedWriter")
        
        
    def __enter__(self):
        return self

    
    def close(self):
        self._file.close()


    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


    def insert_3_bytes(self, barray, i1, i2, idx):
        barray[idx+2] = i2%256 # third byte
        barray[idx+1] = (i2>>8) # first half (from right) of second byte

        barray[idx] = (i1>>4) # first byte
        barray[idx+1] += ((i1%16)<<4) # second half (from right) of second byte


    def insert_1_byte(self, barray, i1, i2, idx):
        barray[idx] = i2
        barray[idx] += (i1<<4)


    def na_to_bytes(self, na: NucleicAcid, name: str, with_struct: bool):
        namelen, slen, plen = len(name), len(na), len(na.pairs)
        
        Nbytes = 2 + 3 + namelen + math.ceil(0.5*slen)
        if with_struct: Nbytes += math.ceil(1.5*plen)
        na_bytes = bytearray([0 for _ in range(Nbytes)])

        # meta block
        na_bytes[1] = Nbytes%256
        na_bytes[0] = Nbytes>>8
        pointer = 2

        self.insert_3_bytes(na_bytes, namelen, slen, pointer)
        pointer += 3

        # name
        for c in name:
            na_bytes[pointer] = ord(c)
            pointer += 1

        # seq
        paired_nbs = set([o for o, _ in na.pairs])
        for i in range(0, slen, 2):
            nb1 = NB_DICT[na.seq[i]] | 8*int(i in paired_nbs and with_struct)
            
            nb2 = 0
            if i<(slen-1):
                nb2 += NB_DICT[na.seq[i+1]] | 8*int((i+1) in paired_nbs and with_struct)
            
            self.insert_1_byte(na_bytes, nb1, nb2, pointer)
            pointer += 1

        # pairs
        if with_struct:
            pairs = na.pairs
            for i in range(0, 2*(plen//2), 2): # iterate even part of sequence
                nb1, nb2 = pairs[i][1], pairs[i+1][1]
                self.insert_3_bytes(na_bytes, nb1, nb2, pointer)
                pointer += 3

            if plen%2!=0:
                nb1 = pairs[-1][1]
                na_bytes[pointer] = nb1>>4
                na_bytes[pointer+1] = ((nb1%16)<<4)

        return na_bytes


    def write(self, na: NucleicAcid, 
              write_struct: bool = True, 
              write_meta: bool = True
              ):

        if len(na)>4095:
            raise ValueError(f"Too long sequence ({len(na)} nb), maximum 4095 nbs supported")
        
        name = na.name
        if write_meta and na.meta:
            str_meta = str(na.meta).replace(' ', '').replace("'", '').strip("{}")
            name = f"{name}{META_SEPARATOR}{str_meta}"

        if len(name)>4095:
            _with_meta = '(with meta json) ' if write_meta and na.meta else ''
            raise ValueError(f"Name length {_with_meta}is too long, maximum 4095 characters supported")
        
        if len(rem:=(set(na.seq) - SUPPORTED_NB_SYMBOLS))!=0:
            raise InvalidSequence(f"Only supported symbols - (A G C U T I N), got {', '.join(tuple(rem))}")
        
        na_bytes = self.na_to_bytes(na, name, write_struct)
        na_bytes = bytes(na_bytes)
        self._file.write(na_bytes)


class bnaRead:

    def __init__(self, file: Union[str, Path, BufferedWriter, BufferedRandom]):
        if isinstance(file, (str, Path)):
            self._file = open(file, 'rb')
        elif isinstance(file, (BufferedWriter, BufferedRandom)):
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

        size_block = int.from_bytes(self._file.read(2), 'big', signed=False)
        if size_block==0: return 0
        na_bytes = self._file.read(size_block-2)

        while na_bytes:
            count += 1
            size_block = int.from_bytes(self._file.read(2), 'big', signed=False)
            if size_block==0: break
            na_bytes = self._file.read(size_block-2)
                
        self._file.seek(tell)
        return count
        
        
    def _iterate(self):
        size_block = int.from_bytes(self._file.read(2), 'big', signed=False)
        if size_block==0: return
        na_bytes = self._file.read(size_block-2)

        while True:
            yield self._make_na(na_bytes)

            size_block = int.from_bytes(self._file.read(2), 'big', signed=False)
            if size_block==0: break
            na_bytes = self._file.read(size_block-2)

        
    def __iter__(self) -> Iterator[NucleicAcid]:
        return self._iterator
                    
                    
    def __next__(self) -> NucleicAcid:
        return next(self._iterator)
    
    
    def read_3_bytes(self, barray, idx):
        i1 = (barray[idx]<<4) + (barray[idx+1]>>4)
        i2 = barray[idx+1]&15
        i2 = (i2<<8) + barray[idx+2]
        return i1, i2

    
    def read_1_byte(self, barray, idx):
        return (barray[idx]>>4), (barray[idx]&15)
    

    def parse_name(self, name_str):
        if META_SEPARATOR not in name_str:
            return name_str, None
        
        name, meta_str = name_str.split(META_SEPARATOR)
        meta = {}
        for kv in meta_str.split(','):
            k, v = kv.split(':')
            meta[k] = v

        return name, meta
    

    def _make_na(self, na_bytes: bytes) -> NucleicAcid:
        pointer = 0
        namelen, slen = self.read_3_bytes(na_bytes, pointer)
        pointer += 3

        # name
        if namelen:
            name_str = []
            for _ in range(namelen):
                name_str.append(chr(na_bytes[pointer]))
                pointer += 1
            name_str = ''.join(name_str)
            name, meta = self.parse_name(name_str)
        else:
            name, meta = None, None

        # seq
        seq = []
        paired_nbs = []
        for i in range(0, 2*(slen//2), 2):
            i1, i2 = self.read_1_byte(na_bytes, pointer)
            
            seq.append(INV_NB_DICT[(i1&7)])
            seq.append(INV_NB_DICT[(i2&7)])

            if (i1&8): paired_nbs.append(i)
            if (i2&8): paired_nbs.append(i+1)

            pointer += 1

        if slen%2!=0:
            i1 = na_bytes[pointer]>>4
            seq.append(INV_NB_DICT[(i1&7)])
            pointer += 1

        seq = ''.join(seq)

        # pairs
        compl_nbs = []
        plen = len(paired_nbs)
        for i in range(0, 2*(plen//2), 2):
            i1, i2 = self.read_3_bytes(na_bytes, pointer)

            compl_nbs.append(i1)
            compl_nbs.append(i2)
            pointer+=3

        if plen%2!=0:
            i1 = (na_bytes[pointer+1]>>4) + (na_bytes[pointer]<<4)
            compl_nbs.append(i1)

        pairs = tuple(zip(paired_nbs, compl_nbs))

        # na
        na = NucleicAcid()
        if name: na.name = name
        if meta: na.meta.update(meta)
        
        for i, nb in enumerate(seq):
            _ = na._add_node(nb)
                
        if len(pairs):
            for o, e in pairs:
                na._add_bond(o, e)
        else:
            na.__dict__['struct'] = None

        return na
    

bnaWrite.__doc__ = format_doc
bnaRead.__doc__ = format_doc