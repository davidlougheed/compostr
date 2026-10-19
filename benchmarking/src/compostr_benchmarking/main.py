import time

from compostr import decomposition
from trviz.decomposer import Decomposer

from .generator import GeneratedSequence, Generator

LENGTHS = [100, 1000, 10000]
SAMPLES = 100


def assess_compostr(seq: GeneratedSequence, motifs: list[str]):
    d = decomposition.MotifSequenceDecomposer([s.encode("ascii") for s in motifs], 4, -5, 7, -1)
    start = time.perf_counter()
    res = d.decompose(bytes(seq))
    duration = time.perf_counter() - start
    return duration


def assess_trviz(seq: GeneratedSequence, motifs: list[str]):
    d = Decomposer()
    start = time.perf_counter()
    res = d.decompose(str(seq), motifs)
    duration = time.perf_counter() - start
    return duration


ASSESS_FUNCTIONS = {
    "compostr": assess_compostr,
    "trviz": assess_trviz,
}

CALLERS = ("compostr", "trviz")


def test_no_interruption():
    records = []
    minor_freq = 0.1

    for length in LENGTHS:
        g = Generator({"CAG": 1.0 - minor_freq, "CCG": minor_freq}, 0.0, (0, 0.0))
        for _ in range(SAMPLES):
            seq = g.generate_sequence(length)
            for caller in CALLERS:
                res = ASSESS_FUNCTIONS[caller](seq, ["CAG", "CCG"])
                records.append({"length": length, "minor_freq": minor_freq, "duration": res, "caller": caller})


def main() -> None:
    # Test cases:

    # Pure tandem repeats of different lengths (mostly a performance test / baseline), possibly split by different
    # lengths and/or motif sizes (STR vs VNTR)

    # Tandem repeats with different levels of alternate motif: 10%, 25%, 50%
    test_no_interruption()

    # Tandem repeats with different levels of interruption: 5%, 10%, ... (although if TRviz doesn't support at all...)

    # Some kind of test of novel motifs...


if __name__ == "__main__":
    main()
