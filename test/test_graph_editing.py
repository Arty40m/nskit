import pytest
from nskit import NA
from nskit.exceptions import InvalidSequence, InvalidStructure



class TestBondEditing:

    @pytest.fixture(scope="session")
    def na(self):
        return NA('..(((..[[...).))..]].') # len - 21
    
    
    @pytest.mark.parametrize(
        "i, compl",
        [
            (0, None), 
            (2, 15), 
            (19, 7),
            (-1, None), 
            (-2, 7)
         ]
    )
    def test_complementary_search(self, i, compl, na):
        assert na.complnb(i)==compl
        
    
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
            # self-bond
            (1, 1),

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

            
    @pytest.mark.parametrize(
        "min_pin_size, target",
        [#       ..((((.))))..((()))..(.)...((.[[.))...]]..
            (1, '..((((.))))..((..))..(.)...((.[[.))...]]..'), 
            (2, '..(((...)))..((..))........((.[[.))...]]..'), 
            (3, '..(((...)))..(....)........((.[[.))...]]..'), 
            (4, '..((.....))..(....)........((.[[.))...]]..'), 
            (5, '..((.....))................(..[[..)...]]..'), 
            (6, '..(.......)................(..[[..)...]]..'), 
            (7, '..(.......)...................(........)..'), 
            (8, '..............................(........)..'), 
            (9, '..........................................'), 
            (10, '..........................................'), 
         ]
    )
    def test_sharp_hairpin_fix(self, min_pin_size, target):
        na = NA('..((((.))))..((()))..(.)...((.[[.))...]]..')
        na.fix_sharp_hairpins(min_pin_size)
        assert na.struct==target
    
    
    
    
    
    
    