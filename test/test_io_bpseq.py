import pytest
import tempfile
from nskit import NA, bpseqRead, bpseqWrite
from nskit.exceptions import InvalidStructure


self_bound_bpseq = \
"""1 G 3
2 G 0
3 G 1
4 G 0
5 G 5
6 G 0
"""

one_direction_bpseq = \
"""1 G 0
2 G 4
3 G 0
4 G 0
"""

two_bonds_bpseq = \
"""1 G 3
2 G 0
3 G 5
4 G 0
5 G 0
"""

out_of_range_bond_bpseq = \
"""1 G 0
2 G 0
3 G 6
4 G 0
"""

negative_bond_bpseq = \
"""1 G 0
2 G 0
3 G -2
4 G 0
"""

unconsistent_bpseq = \
"""1 G 0
2 G 0
5 G 0
6 G 0
"""


class TestBpseq:

    @pytest.mark.parametrize(
        "na",
        [
            NA('AAA', '...', name='Seq1'), 
            NA('UUUUCCCC', '((...)).', name='Seq2'),
            NA('CGCGCGCGCGCGCGCGCGCGCGCAG', '..(((..(..)..[..)))...]..', name='Seq3', meta={'param1':'1', 'param2':'2'}),
            NA('ATATCGCGATAC', '(.[.{.)]..}.', name='Seq4', meta={'param':'[1,2,3]'})
         ]
    )
    def test_read_write(self, na):
        fp = tempfile.TemporaryFile('w+')
        with bpseqWrite(fp) as w:
            w.write(na)

            fp.seek(0)

            with bpseqRead(fp) as f:
                na_ = f.read()

            assert na.seq==na_.seq
            assert na.struct==na_.struct
            assert na.meta==na_.meta

    
    def test_self_bound(self):
        fp = tempfile.TemporaryFile('w+')
        fp.write(self_bound_bpseq)
        fp.seek(0)

        with bpseqRead(fp, raise_na_errors=True) as f:
            with pytest.raises(InvalidStructure):
                _ = f.read()


    def test_one_direction(self):
        fp = tempfile.TemporaryFile('w+')
        fp.write(one_direction_bpseq)
        fp.seek(0)

        with bpseqRead(fp, raise_na_errors=True) as f:
            with pytest.raises(InvalidStructure):
                _ = f.read()

    
    def test_two_bonds(self):
        fp = tempfile.TemporaryFile('w+')
        fp.write(two_bonds_bpseq)
        fp.seek(0)

        with bpseqRead(fp, raise_na_errors=True) as f:
            with pytest.raises(InvalidStructure):
                _ = f.read()


    def test_out_of_range_bond(self):
        fp = tempfile.TemporaryFile('w+')
        fp.write(out_of_range_bond_bpseq)
        fp.seek(0)

        with bpseqRead(fp, raise_na_errors=True) as f:
            with pytest.raises(InvalidStructure):
                _ = f.read()


    def test_negative_bond(self):
        fp = tempfile.TemporaryFile('w+')
        fp.write(negative_bond_bpseq)
        fp.seek(0)

        with bpseqRead(fp, raise_na_errors=True) as f:
            with pytest.raises(InvalidStructure):
                _ = f.read()

                
    def test_unconsistent(self):
        fp = tempfile.TemporaryFile('w+')
        fp.write(unconsistent_bpseq)
        fp.seek(0)

        with bpseqRead(fp, raise_na_errors=True) as f:
            with pytest.raises(InvalidStructure):
                _ = f.read()


    
    
    
    
    