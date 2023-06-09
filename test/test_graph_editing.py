import pytest
from nskit import NA
from nskit.exceptions import InvalidSequence, InvalidStructure



class TestBondEditing:

    @pytest.mark.parametrize(
        "struct, i, j, target",
        [
            ('.((..)).....', 8, 11, '.((..)).(..)'), 
            ('.((..)).....', 4, 8, '.((.[)).]...'), 
            ('.((..)).....', -4, 11, '.((..)).(..)'), 
         ]
    )
    def test_bond_addition_deletion(self, struct, i, j, target):
        na = NA(struct)
        na.join(i, j)
        assert na.struct==target

        na.split(i, j)
        assert na.struct==struct


    @pytest.fixture(scope="session")
    def na(self):
        return NA('..(((..[[...).))..]].') # len - 21
    

    @pytest.mark.parametrize(
        "i, j",
        [
            (0, 22), 
            (-1, -22)
         ]
    )
    def test_join_out_of_range(self, i, j, na):
        with pytest.raises(IndexError):
            na.join(i, j)


    @pytest.mark.parametrize(
        "i, j",
        [
            # sharp
            (1, 1),
            (0, 1),

            # multibond
            (0, 2),
            (2, 6), 
         ]
    )
    def test_join_value_error(self, i, j, na):
        with pytest.raises(ValueError):
            na.join(i, j)


    @pytest.mark.parametrize(
        "i, j",
        [
            (0, 0), 
            (2, 2),
            (0, 7),
            (7, 8),
            (0, -1),
            
         ]
    )
    def test_split_no_bond(self, i, j, na):
        with pytest.raises(ValueError):
            na.split(i, j)

    
    