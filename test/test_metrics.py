import pytest
import nskit as nsk
import numpy as np



class TestBinaryMetrics:

    def test_length_error(self):
        with pytest.raises(ValueError):
            _ = nsk.metrics.binary_eval(".((..)).", "((..))")
            
            
    def test_empty_error(self):
        with pytest.raises(ValueError):
            _ = nsk.metrics.binary_eval("......", "((..))")
        
        
    @pytest.mark.parametrize(
        "true, pred, y",
        [
            ("((((...))))", "((((...))))", 1.), 
            ("((((...))))", "(((((.)))))", 1.), 
            ("((((...))))", ".(((...))).", 0.75), 
            ("..((..))..", "..((....))", 0.), 
            ("((.[[.))..]]", "((....))....", 0.5), 
            ("..((...))..", "...........", 0.)
         ]
    )
    def test_recall(self, true, pred, y):
        x = nsk.metrics.recall(true, pred)
        assert np.isclose(x, y)
        
    
    @pytest.mark.parametrize(
        "true, pred, y",
        [
            ("((..))", "((..))", 1.), 
            ("(((..)))", "((....))", 1.), 
            ("(((..)))", "(((.).))", 2/3), 
            ("..((...))..", "((.....))..", 0.), 
            ("..((...))..", "...........", 0.)
         ]
    )
    def test_precision(self, true, pred, y):
        x = nsk.metrics.precision(true, pred)
        assert np.isclose(x, y)