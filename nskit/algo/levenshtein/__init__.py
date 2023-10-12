from ._levenshtein import c_levenshtein
from typing import Union
from ...containers.nucleic_acid import NucleicAcid



def levdist(a: Union[str, NucleicAcid], 
            b: Union[str, NucleicAcid], 
            ins: float = 1., 
            rm: float = 1., 
            sub: float = 1.
           ) -> float:
    """
    Calculates levenshtein distance between two strings or NucleicAcid sequences.

    :param a: first ascii string or NucleicAcid.
    :param b: second ascii string or NucleicAcid.
    :param ins: insert weight.
    :param rm: delete(remove) weight.
    :param sub: substitute weight.

    :return: distance float value.
    """
    
    if isinstance(a, NucleicAcid):
        a = a.seq
    if isinstance(b, NucleicAcid):
        b = b.seq
        
    a = a.encode('ascii')
    b = b.encode('ascii')
    
    if len(a)==0 or len(b)==0:
        return float((len(a) | len(b))*ins)

    dist = c_levenshtein(a, b, 
                         len(a), len(b), 
                         float(ins), 
                         float(rm), 
                         float(sub)
                        )

    return dist


__all__ = ["levdist"]