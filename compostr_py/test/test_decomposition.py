from compostr import decomposition


def test_init_decomposer():
    decomposition.MotifSequenceDecomposer([b"CAG", b"CCG"], 4, -5, 7, -1)


def test_basic_decomposition():
    d = decomposition.MotifSequenceDecomposer([b"CAG", b"CCG"], 4, -5, 7, -1)
    seq = b"CAGCAGCCGCAG"
    res = d.decompose(seq)
    assert isinstance(res, decomposition.MotifSequenceDecomposition)
    assert res.sequence_items(seq, False) == [b"CAG", b"CAG", b"CCG", b"CAG"]
