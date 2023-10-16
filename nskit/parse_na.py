from collections import defaultdict
from typing import Optional, Union, Tuple, List

from .containers import NucleicAcid
from .exceptions import InvalidSequence, InvalidStructure



DOT_STRUCTURE_SYMBOLS = set(".()[]{}<>AaBbCcDdEeFf")
STRUCTURE_DETECTION_SYMBOLS = set(".()[]{}<>")


class NaStack:
    inv_dict = {')':'(', ']':'[', '}':'{', '>':'<', 'a':'A', 
                'b':'B', 'c':'C', 'd':'D', 'e':'E', 'f':'F'}
    
    def __init__(self):
        self.st = defaultdict(list)

    def __setitem__(self, k, v):
        self.st[k].append(v)
    
    def __getitem__(self, k):
        if len(self.st[self.inv_dict[k]]):
            return self.st[self.inv_dict[k]].pop()
        return None
        
    def isempty(self):
        for k in self.st:
            if len(self.st[k])!=0:
                return False
        return True


def parse_arguments(a:Optional[str], b:Optional[str]) -> Tuple[Optional[str], Optional[str]]:
    if b is not None:
        seq = a
        struct = b
        
        if len(seq)!=len(struct):
            raise InvalidStructure(f"Sequence and dot structure must have equal length, "
                                   f"got {len(seq)} and {len(struct)}")
        
    elif len(set(a) & STRUCTURE_DETECTION_SYMBOLS):
        struct = a
        seq = None
    else:
        seq = a
        struct = None
        
    return seq, struct


def parse_structure(struct: str, 
                    ignore_unclosed_bonds: bool, 
                   ) -> List:
    pairs = []
    if len(rem:=(set(struct) - DOT_STRUCTURE_SYMBOLS))!=0:
        raise InvalidStructure(f"Dot structure contains invalid symbols - {', '.join(tuple(rem))}")

    stack = NaStack()
    for i, s in enumerate(struct):
        if s == '.':
            continue

        if s in '([{<ABCDEF':
            stack[s] = i

        elif s in ')]}>abcdef':
            op_idx = stack[s]

            if op_idx is None:
                if ignore_unclosed_bonds:
                    continue
                raise InvalidStructure(f"Closing bond {s} at index {i} has no open pair, "
                                       f"use ignore_unclosed_bonds=True to omit such bonds")
                    
            pairs.append((op_idx, i))

    if not stack.isempty() and not ignore_unclosed_bonds:
        raise InvalidStructure(f"Structure contains unclosed bonds, use ignore_unclosed_bonds=True to omit such bonds")
        
    return pairs
    
    
def NA(a: Union[str, NucleicAcid], b: Optional[str] = None, /, *, 
       name: Optional[str] = None, 
       meta: Optional[dict] = None,
       ignore_unclosed_bonds: bool = False, 
       upper_sequence: bool = False,
      ) -> NucleicAcid:
    """
    Parse dotbracket strings into NucleicAcid.

    :param a: sequence or structure if sequence is not provided.
    :param b: structure if sequence is provided.
    :param name: na name.
    :param meta: dictionary of meta information convertable to string.
    :param ignore_unclosed_bonds: omit single unpaired parentheses without raising error. Default - False.
    :param upper_sequence: upper sequence characters. Default - False.

    :return: NucleicAcid object.
    """ 
    
    if not a: 
        raise ValueError("Empty data")
        
    if isinstance(a, NucleicAcid): 
        return a
        
    seq, struct = parse_arguments(a, b)
     
    # validate sequence
    if seq: 
        if not seq.isalpha():
            raise InvalidSequence(f"Sequence must contain only alphabetic characters")
            
        if upper_sequence: 
            seq = seq.upper()
    else:
        seq = 'N'*len(struct)
        
    # parse structure
    if struct:
        pairs = parse_structure(struct, ignore_unclosed_bonds)
    
    # create graph
    na = NucleicAcid()
    if name: na.name = name
    if meta: na.meta.update(meta)
    
    for i, nb in enumerate(seq):
        _ = na._add_node(nb)
            
    if struct:
        for o, e in pairs:
            na._add_bond(o, e)
    else:
        na.__dict__['struct'] = None
    
    return na