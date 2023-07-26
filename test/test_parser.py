import pytest
from nskit import NA
from nskit.exceptions import InvalidSequence, InvalidStructure



class TestSequence:

    @pytest.mark.parametrize(
        "data",
        ['', None]
    )
    def test_empty_sequnce(self, data):
        with pytest.raises(ValueError):
            _ = NA(data)
        

    @pytest.mark.parametrize(
        "data",
        ['56', '%AAA', '-UU']
    )
    def test_alphabetic(self, data):
        with pytest.raises(InvalidSequence):
            _ = NA(data)


class TestStructure:

    @pytest.mark.parametrize(
        "data",
        ['R..r', '.(.--.).', '.AAA']
    )
    def test_unknown_symbols(self, data):
        with pytest.raises(InvalidStructure):
            _ = NA(data)
            
            
    @pytest.mark.parametrize(
        "inp, out",
        [
            ('.[[..]]..', '.((..))..'), 
            ('..((..).', '...(..).'),
            ('..(..)).', '..(..)..'),
            ('.((.<.]}.>', '....(....)'),
            
        ]
    )
    def test_unclosed_bonds_process(self, inp, out):
        na = NA(inp, ignore_unclosed_bonds=True)
        assert na.struct==out
        
        
    @pytest.mark.parametrize(
        "struct",
        ['....', '..((..]].', '..).']
    )
    def test_empty_structure(self, struct):
        na = NA(struct, ignore_unclosed_bonds=True)
        assert na.struct == '.'*len(na)


    @pytest.mark.parametrize(
        "struct",
        ['((()))', '..((..)).()', '..{{}}.']
    )
    def test_sharp_helix_error(self, struct):
        with pytest.raises(InvalidStructure):
            _ = NA(struct)


    @pytest.mark.parametrize(
        "struct, fixed",
        [
            ('((()))', '((..))'), 
            ('..().().', '........'), 
            ('.(.).[[]].', '.(.).(..).'), 
        ]
    )
    def test_sharp_helix_fix(self, struct, fixed):
        na = NA(struct, fix_sharp_helixes=True)
        assert na.struct==fixed
        
    
    @pytest.mark.parametrize(
        "struct, fixed",
        [
            ('((()))', '((()))'), 
            ('..((..)).()', '..((..)).()'), 
            ('..{{}}.', '..(()).')
        ]
    )
    def test_allow_sharp_helix(self, struct, fixed):
        na = NA(struct, allow_sharp_helixes=True)
        assert na.struct==fixed
        

class TestGraphParse:

    @pytest.mark.parametrize(
        "struct, is_knot",
        [
            ("..((..[[.))..]]..", True), 
            ("...((..))..[[[..]]]", False),
            ("((..(((..)))...((...))..))", False), 
            ("....", False)
         ]
    )
    def test_is_knot(self, struct, is_knot):
        assert NA(struct).is_knot()==is_knot


    @pytest.mark.parametrize(
        "struct, orders",
        [
            ("..((..[[.))..]]..", (0, 1)), 
            ("...((..))..[[[..]]]", (0, 0)),
            ("..(((.[[.)))..{{{.((..]].))..}}}.", (0, 1, 0, 0)),
            (".((.[[.)).]]..([{<)]}>....", (0, 1, 0, 1, 2, 3)) 
         ]
    )
    def test_helix_orders(self, struct, orders):
        na = NA(struct)
        for i, order in enumerate(na.helix_orders):
            assert order == orders[i]
    

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        