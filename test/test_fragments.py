import pytest
from nskit import NA
from nskit.containers import Helix, Hairpin, InternalLoop, Bulge, Junction

class TestLoops:
    
    @pytest.mark.parametrize(
        "na, target_fragments",
        [
            (NA('.((.((..((((..))..))........((..((.[[.))))....]]...))..))...'), 
             (InternalLoop, Junction, Bulge, Hairpin, Bulge, Hairpin)), 
            
            (NA('..((..(((((....(((..((..[[[..))..))))))...))..))..]]]..'), 
             (InternalLoop, Bulge, Bulge, InternalLoop, Hairpin)), 
         ]
    )
    def test_fragments_order(self, na, target_fragments):
        for i, l in enumerate(na.loops):
            assert isinstance(l, target_fragments[i])
            
    
    def test_sharp_hairpin(self):
        na = NA('..()..')
        assert len(na.loops[0].nts)==2
        
        
    def test_no_loops(self):
        na = NA('.....')
        assert len(na.loops)==0
        
        
    def test_hairpin_length(self):
        na = NA('..(((.....)))..')
        assert len(na.hairpins[0].nts)==7
        assert len(na.hairpins[0])==7
            
            
    def test_internal_loop_length(self):
        na = NA('..(((...((...))...)))..')
        assert len(na.internal_loops[0].nts)==10
        assert len(na.internal_loops[0])==10
            
            
    def test_bulge_length(self):
        na = NA('..(((((...))..)))..')
        assert len(na.bulges[0].nts)==6
        assert len(na.bulges[0])==6
            
            
    def test_junction_length(self):
        na = NA('..(((..((...))...((...))...((..))...)))..')
        assert len(na.junctions[0].nts)==19
        assert len(na.junctions[0])==19
            
          
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
        
        
class TestKnotsInLoop:
    
    def test_not_knots(self):
        na = NA('.((..))..')
        assert na.loops[0].has_knot() == False
        
        
    @pytest.mark.parametrize(
        "na, knots",
        [
            (NA('.((.[[.)).]].'), ((4, 5),)), 
            (NA('.((.[[.[[.)).]].]].'), ((4, 5), (7, 8))), 
            (NA('.((.[[.{{.)).]].}}.'), ((4, 5), (7, 8))), 
            (NA('.((.[[{{..)).]].}}.'), ((4, 5), (6, 7))), 
         ]
    )
    def test_hairpin_knots(self, na, knots):
        assert na.loops[0].has_knot()
        assert na.loops[0].knots == knots
    
    
    def test_junction_knots(self):
        na = NA('..((..((...))..((..[[.)).]].))..')
        assert na.junctions[0].knots == ((25, 26),)
        assert na.hairpins[1].knots == ((19, 20),)
    
    
class TestDanglingEnds:
    
    @pytest.mark.parametrize(
        "na, ends",
        [
            (NA('.((...))..'), ((0,), (8, 9)) ), 
            (NA('((...))..'), (tuple(), (7, 8)) ), 
            (NA('..((..))'), ((0, 1), tuple()) ), 
            (NA('((..))'), (tuple(), tuple()) ), 
            (NA('....'), ((0, 1, 2, 3), tuple()) ), 
         ]
    )
    def test_ends(self, na, ends):
        assert na.dangling_ends==ends
    
    
class TestfragmentProperties:
    
    def test_symmetric_internal_loop(self):
        na = NA('((..((...))..))')
        assert na.internal_loops[0].is_symmetric()
        assert na.internal_loops[0].short_side==2
        
    
    @pytest.mark.parametrize(
        "na, short_side, long_side",
        [
            (NA('((.((...))..))'), 1, 2), 
            (NA('((....((...))..))'), 2, 4), 
            (NA('((..((..[[.))..]]...))'), 2, 7), 
         ]
    )
    def test_asymmetric_internal_loop(self, na, short_side, long_side):
        assert not na.internal_loops[0].is_symmetric()
        assert na.internal_loops[0].short_side==short_side
        assert na.internal_loops[0].long_side==long_side
        
        
        
    
    
        
        
        