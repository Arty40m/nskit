from typing import Union
from pathlib import Path
from io import BufferedWriter
import math

from ..nucleic_acid import NucleicAcid
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


class bnaWrite():
    
    def __init__(self, file: Union[str, Path, BufferedWriter], *, 
                 append: bool = False,
                ):
        
        if isinstance(file, (str, Path)):
            self._file = open(file, 'ab' if append else 'wb')
        elif isinstance(file, BufferedWriter):
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
            str_meta = str(na.meta).replace(' ', '').replace("'", '')
            name = f"{name}{META_SEPARATOR}{str_meta}"

        if len(name)>4095:
            _with_meta = '(with meta json) ' if write_meta and na.meta else ''
            raise ValueError(f"Name length {_with_meta}is too long, maximum 4095 characters supported")
        
        if len(rem:=(set(na.seq) - SUPPORTED_NB_SYMBOLS))!=0:
            raise InvalidSequence(f"Only supported symbols - (A G C U T I N), got {', '.join(tuple(rem))}")
        
        na_bytes = self.na_to_bytes(na, name, write_struct)
        na_bytes = bytes(na_bytes)
        self._file.write(na_bytes)


bnaWrite.__doc__ = format_doc