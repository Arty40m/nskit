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
        "data",
        ['....', '..((..]].', '..).']
    )
    def test_empty_structure(self, data):
        with pytest.raises(InvalidStructure):
            _ = NA(data, filter_linear_structures=True, ignore_unclosed_bonds=True)
            
    

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        