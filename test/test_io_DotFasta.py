import pytest
import tempfile
from nskit import NA, DotRead, DotWrite, FastaRead, FastaWrite
from nskit.exceptions import InvalidFasta



dot_file = \
""">Seq1

AAA
>Seq2
UUUU
((.)

>Seq3
CCCGGG
param1: 1
param2: 2
>Seq4
AAATTT
(.[.)]
param: [1,2,3]
"""

nas = [
    NA('AAA', name='Seq1'), 
    NA('UUUU', '.(.)', name='Seq2'),
    NA('CCCGGG', name='Seq3', meta={'param1':'1', 'param2':'2'}),
    NA('AAATTT', '(.[.)]', name='Seq4', meta={'param':'[1,2,3]'})
]

def file_with_sharp_helix():
    fp = tempfile.TemporaryFile('w+')
    fp.write(">seq\nCUUCGUGGC\n.().((.))")
    fp.seek(0)
    return fp


class TestDotBracket:

    def test_read(self):
        fp = tempfile.TemporaryFile('w+')
        fp.write(dot_file)
        fp.seek(0)

        with DotRead(fp, ignore_unclosed_bonds=True) as f:
            for i, na in enumerate(f):
                true_na = nas[i]
                assert na.seq==true_na.seq
                assert na.struct==true_na.struct
                assert na.name==true_na.name
                assert na.meta==true_na.meta

    
    def test_sharp_helix(self):
        fp = file_with_sharp_helix()
        with DotRead(fp) as f:
            na = next(f)
            assert na is None
            
            
    def test_sharp_helix_fix(self):
        fp = file_with_sharp_helix()
        with DotRead(fp, fix_sharp_helixes=True) as f:
            na = next(f)
            assert na.struct == "....((.))"
            
            
    def test_allow_sharp_helix(self):
        fp = file_with_sharp_helix()
        with DotRead(fp, allow_sharp_helixes=True) as f:
            na = next(f)
            assert na.struct == ".().((.))"
                
                
    def test_empty_na(self):
        fp = tempfile.TemporaryFile('w+')
        fp.write(">seq1\n\n>seq2\naaa")
        fp.seek(0)

        with DotRead(fp) as f:
            with pytest.raises(InvalidFasta):
                _=next(f)


    def test_write(self):
        fp = tempfile.TemporaryFile('w+')
        with DotWrite(fp) as w:
            for na in nas:
                w.write(na)

            fp.seek(0)

            with DotRead(fp) as f:
                for i, na in enumerate(f):
                    true_na = nas[i]
                    assert na.seq==true_na.seq
                    assert na.struct==true_na.struct
                    assert na.name==true_na.name
                    assert na.meta==true_na.meta


fasta_file = \
""">seq
;LCBO - Prolactin precursor - Bovine
UAUGUGAAUCUCCUACGUGAACGUUGUGCCAAUUUUGCCUCGUCGCCGCCCCCCGAAGGCCAUGACCUGUUGGUACACUA
CAUUGGUUUUGCACCAUAAUGUUCCCUCCUCGGGUAUAUCACGGUGUUCAUCGCCACAGGAAAUUUCCGUAAGUUCCAUA
GUCGCCGCCCCCCGAAGGCCAUGACCUGUUGGUACACUA*
"""
fasta_seq = "UAUGUGAAUCUCCUACGUGAACGUUGUGCCAAUUUUGCCUCGUCGCCGCCCCCCGAAGGCCAUGACCUGUUGGUACACUA" + \
"CAUUGGUUUUGCACCAUAAUGUUCCCUCCUCGGGUAUAUCACGGUGUUCAUCGCCACAGGAAAUUUCCGUAAGUUCCAUA" + \
"GUCGCCGCCCCCCGAAGGCCAUGACCUGUUGGUACACUA"

class TestFasta:
    
    def test_read(self):
        fp = tempfile.TemporaryFile('w+')
        fp.write(fasta_file)
        fp.seek(0)

        with FastaRead(fp) as f:
            na = next(f)
            assert na.seq==fasta_seq
            assert na.name=="seq"


    def test_write(self):
        na1 = NA(fasta_seq)

        fp = tempfile.TemporaryFile('w+')
        with FastaWrite(fp) as w:
            w.write(['test_seq', fasta_seq])
            
            fp.seek(0)
            assert len(fp.readlines())==4
            fp.seek(0)

            with FastaRead(fp) as f:
                na2 = next(f)
            
        assert na1.seq==na2.seq