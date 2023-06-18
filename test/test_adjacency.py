import pytest
import numpy as np
from nskit import NA, NucleicAcid
from nskit.exceptions import InvalidStructure, InvalidAdjacency


asymmetric_adjacency = np.array([
    [0, 0, 1], 
    [0, 0, 0], 
    [0, 0, 0],  
])

values_bt1 = np.array([
    [0, 0, 2], 
    [0, 0, 0], 
    [2, 0, 0], 
])

values_lt0 = np.array([
    [0, 0, -1], 
    [0, 0, 0], 
    [-1, 0, 0], 
])

values_float = np.array([
    [0, 0, 1.1], 
    [0, 0, 0], 
    [1.1, 0, 0], 
])

multiple_bonds = np.array([
    [0, 0, 1, 0, 1], 
    [0, 0, 0, 0, 0], 
    [1, 0, 0, 0, 0], 
    [0, 0, 0, 0, 0], 
    [1, 0, 0, 0, 0], 
])

sharp_helix = np.array([
    [0, 0, 0, 1], 
    [0, 0, 1, 0], 
    [0, 1, 0, 0], 
    [1, 0, 0, 0], 
])

class TestAdjacency:

    @pytest.mark.parametrize(
        "struct",
        [
            '....', 
            '..((....))..',
            '..(((...))..).',
            '..((.[[.)).].]',
            '([{)]}',
            '(((..)))'
         ]
    )
    def test_recreation(self, struct):
        na1 = NA(struct)
        na2 = NucleicAcid.from_adjacency(na1.get_adjacency())
        assert na1.struct==na2.struct


    def test_asymmetry(self):
        with pytest.raises(InvalidAdjacency):
            _ = NucleicAcid.from_adjacency(asymmetric_adjacency)

    @pytest.mark.parametrize(
        "adj",
        [
            values_bt1,
            values_lt0,
            values_float
         ]
    )
    def test_values(self, adj):
        with pytest.raises(InvalidAdjacency):
            _ = NucleicAcid.from_adjacency(adj)


    def test_multiple_bonds(self):
        with pytest.raises(InvalidAdjacency):
            _ = NucleicAcid.from_adjacency(multiple_bonds)


    def test_sharp_helix_error(self):
        with pytest.raises(InvalidStructure):
            _ = NucleicAcid.from_adjacency(sharp_helix)

    
    def test_sharp_helix_fix(self):
        na = NucleicAcid.from_adjacency(sharp_helix, fix_sharp_helixes=True)
        assert na.struct == '(..)'