import pytest
import tempfile
from nskit import NA, dotRead, dotWrite, fastaRead, fastaWrite, dotLinesRead, dotLinesWrite
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


class TestDotBracket:

    def test_read(self):
        fp = tempfile.TemporaryFile('w+')
        fp.write(dot_file)
        fp.seek(0)

        with dotRead(fp, ignore_unclosed_bonds=True) as f:
            for i, na in enumerate(f):
                true_na = nas[i]
                assert na.seq==true_na.seq
                assert na.struct==true_na.struct
                assert na.name==true_na.name
                assert na.meta==true_na.meta

                
    def test_empty_na(self):
        fp = tempfile.TemporaryFile('w+')
        fp.write(">seq1\n\n>seq2\naaa")
        fp.seek(0)

        with dotRead(fp) as f:
            with pytest.raises(InvalidFasta):
                _=next(f)


    def test_write(self):
        fp = tempfile.TemporaryFile('w+')
        with dotWrite(fp) as w:
            for na in nas:
                w.write(na)

            fp.seek(0)

            with dotRead(fp) as f:
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

        with fastaRead(fp) as f:
            na = next(f)
            assert na.seq==fasta_seq
            assert na.name=="seq"


    def test_write(self):
        na1 = NA(fasta_seq, name='test_seq')

        fp = tempfile.TemporaryFile('w+')
        with fastaWrite(fp) as w:
            w.write(na1)
            
            fp.seek(0)
            assert len(fp.readlines())==4
            fp.seek(0)

            with fastaRead(fp) as f:
                na2 = next(f)
            
        assert na1.seq==na2.seq

        
raw_lines = [
    ('>seq1', 'fff', '...'), 
    ('>seq1', 'ffttt', '...(())', 'e:2'), 
    ('>seq1', 'CCCC')
]

class TestFasta:
    
    def test_read_write(self):
        fp = tempfile.TemporaryFile('w+')
        with dotLinesWrite(fp) as w:
            for lines in raw_lines:
                w.write(lines)

            fp.seek(0)

            with dotLinesRead(fp) as f:
                for i, lines in enumerate(f):
                    assert lines==raw_lines[i]
            
            