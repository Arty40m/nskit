import pytest
from nskit import NA
from nskit.containers import Helix, Hairpin, InternalLoop, Bulge, Junction

class TestLoops:
    
    @pytest.mark.parametrize(
        "na, target_parts",
        [
            (NA('.((.((..((((..))..))........((..((.[[.))))....]]...))..))...'), 
             (InternalLoop, Junction, Bulge, Hairpin, Bulge, Hairpin)), 
            
            (NA('..((..(((((....(((..((..[[[..))..))))))...))..))..]]]..'), 
             (InternalLoop, Bulge, Bulge, InternalLoop, Hairpin)), 
         ]
    )
    def test_parts_order(self, na, target_parts):
        for i, l in enumerate(na.loops):
            assert isinstance(l, target_parts[i])
            
            
    def test_no_loops(self):
        na = NA('.....')
        assert len(na.loops)==0
        
        
    def test_hairpin_length(self):
        na = NA('..(((.....)))..')
        assert len(na.loops[0].loop)==7
        assert len(na.loops[0])==7
            
            
    def test_internal_loop_length(self):
        na = NA('..(((...((...))...)))..')
        assert len(na.loops[0].loop)==10
        assert len(na.loops[0])==10
            
            
    def test_bulge_length(self):
        na = NA('..(((((...))..)))..')
        assert len(na.loops[0].loop)==6
        assert len(na.loops[0])==6
            
            
    def test_junction_length(self):
        na = NA('..(((..((...))...((...))...((..))...)))..')
        assert len(na.loops[0].loop)==19
        assert len(na.loops[0])==19
            
          
    @pytest.mark.parametrize(
        "na, n",
        [
            (NA('..(((..((...))...((...))...((..))...)))..'), 4), 
            (NA('..(((...)))..'), 1), 
            (NA('..(((..((...))..)))..'), 2), 
            (NA('..(((..((...)))))..'), 2), 
         ]
    )
    def test_number_of_branches(self, na, n):
        assert len(na.loops[0].branches)==n
        
        
        
        
        
        