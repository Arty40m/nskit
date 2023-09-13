import pytest
from nskit import NA
from nskit.algo import levdist



class TestLevenshtein:

    @pytest.mark.parametrize(
        "a, b, dist",
        [
            ('hodqdbhoEbao', 'hodqdbhoEbao', 0.),
            ('hodqdbhoEbao', 'E', 11.),
            ('hodqdbhoEbao', '', 12.),
            ('E', 'ER', 1.),
            (
                ('UAGCGCGGAGCUACCUGAUAUGGUAUUACGCGCAAGACGAUCGUUCGCUUCGGUACGGA'
                'AAGUUCAUCGUGCAAACCUGCCGCGAUACAUUGAUUGCAUGUAAACGUUUUCCACUAG'
                'UGCUUCUAUUCCUCGCCACACUUCCCUGUUUCUCCCUCUCAAUAUUUCGAGAGCCUCGUU'
                'UUCUAUGGCACGCCUGUCAUGGCGUUCCCCGAAUGGUGCGCCUAUAGGCCGCCAGGACCUA'
                'GUAAAAACUGAACGUUGUAGACAGAGAUAGUGACGCCCAUCGUUGUCAUGGACAGGAUAACU'), 
                
                ('CUACCUCUCAUGGUCCUCGUAGAGACGGUAACAUAGCAGAGAUCCAUACGG'
                'AACACCGGCCAAUAAAUUAAAAAAAACCUGUAGCUUCGCUGACACCGAACU'
                'CGAACGGGAUCCGUUUGAAAGACGACUCUACAAAAUUAUUGGACGGUAA'
                'UCACACACACGGGCUAUCAUUAAGGUCGCUGCGGGCGUCCAGAGUACGC'
                'UCUCGUUAACACCGCCACCGGGAGGGAAAAUC'), 
                
                156.
            )
        ]
    )
    def test_uniform(self, a, b, dist):
        assert levdist(a, b) == dist
        

    @pytest.mark.parametrize(
        "a, b, insert, delete, substitute, dist",
        [
            ('hodqdbhoEbao', 'hodqdbhoEbao', 1, 1, 1, 0.),
            ('qqq', 'qqyq', 2, 100, 1, 2.),
            ('qqyq', 'qqq', 2, 100, 1, 100.),
            ('qqq', 'qqE', 2, 3, 10, 5.),
            ('qqq', 'qqE', 20, 30, 10, 10.),
            ('qqq0834nlosnduuion', 'eiyNLOBLSKNF', 3, 1, 6, 50.),
            ('qqq0834nlosnduuion', 'eiyNLOBLSKNF', .3, 1.5, 1.3, 24.6),
        ]
    )
    def test_weighted(self, a, b, insert, delete, substitute, dist):
        assert levdist(a, b, insert, delete, substitute) == dist
        
    
    @pytest.mark.parametrize(
        "a, b",
        [
            ("AUUGC", "GCGU"),
            (NA("AA"), NA("CCGGC")),
            (NA("CCU"), "ACG"),
            ("CCC", NA("UUU"))
        ]
    )
    def test_NA_str(self, a, b):
        _ = levdist(a, b)
        
        
    @pytest.mark.parametrize(
        "a, Error",
        [
            (None, AttributeError), 
            (chr(1000), UnicodeEncodeError), 
            (chr(200), UnicodeEncodeError), 
            (1, AttributeError), 
            (2., AttributeError), 
            ((2,3), AttributeError)
        ]
    )
    def test_error_data(self, a, Error):
        with pytest.raises(Error):
            _ = levdist(a, "ACGU")


        
    

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        