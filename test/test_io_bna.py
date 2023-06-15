import pytest
import tempfile
from nskit import NA, bnaRead, bnaWrite



nas = [
    NA('AAA', name='Seq1'), 
    NA('UUUU', '.(.)', name='Seq2'),
    NA('CCCGGG', name='Seq3', meta={'param1':'1', 'param2':'2'}),
    NA('AAATTTU', '(.[.)].', meta={'param1':'1', 'param2':'44'})
]


class TestBna:

    def test_read_write(self):
        fp = tempfile.TemporaryFile('w+b')
        with bnaWrite(fp) as w:
            for na in nas:
                w.write(na)
            
            fp.seek(0)

            bnas = []
            with bnaRead(fp) as f:
                _ = len(f)
                for bna in f:
                    _ = len(f)
                    bnas.append(bna)

        for na, bna in zip(nas, bnas):
            assert na.seq == bna.seq
            assert na.struct == bna.struct
            assert na.name == bna.name
            assert na.meta == bna.meta


    def test_write_only_seq(self):
        na = NA('AAATTTU', '(.[.)].', meta={'param1':'1', 'param2':'44'})

        fp = tempfile.TemporaryFile('w+b')
        with bnaWrite(fp) as w:
            w.write(na, write_struct=False, write_meta=False)
            
            fp.seek(0)

            with bnaRead(fp) as f:
                bna = next(f)

            assert na.seq == bna.seq
            assert bna.struct is None
