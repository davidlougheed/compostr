import time

import numpy as np
import polars as pl
from compostr import decomposition
from tqdm import tqdm
from trviz.decomposer import Decomposer

from .generator import GeneratedSequence, Generator, random_dna

LENGTHS = [100, 1000, 10000]
SAMPLES = 100


def assess_compostr(seq: GeneratedSequence, motifs: list[str]):
    d = decomposition.MotifSequenceDecomposer([s.encode("ascii") for s in motifs], 4, -5, 7, False, -1)
    seq_bytes = bytes(seq)
    start = time.perf_counter()
    res = d.decompose(seq_bytes)
    res_seq = res.sequence_items(seq_bytes, False)
    # print("compostr", res_seq)
    duration = time.perf_counter() - start
    return duration


def assess_trviz(seq: GeneratedSequence, motifs: list[str]):
    d = Decomposer()
    seq_str = str(seq)
    start = time.perf_counter()
    res = d.decompose(seq_str, motifs)
    # print("trviz", res)
    duration = time.perf_counter() - start
    return duration


ASSESS_FUNCTIONS = {
    "compostr": assess_compostr,
    "trviz": assess_trviz,
}

CALLERS = ("compostr", "trviz")


def generate_motif_df(motif_freqs: dict[str, float], minor_freq: float, repeat_type: str):
    records = []

    motifs = list(motif_freqs.keys())

    for length in LENGTHS:
        g = Generator(motif_freqs, 0.00, (10, 0.5))
        for _ in tqdm(range(SAMPLES), desc=f"{repeat_type} | no_interruption; length: {length}"):
            seq = g.generate_sequence(length)
            for caller in CALLERS:
                res = ASSESS_FUNCTIONS[caller](seq, motifs)
                records.append({"length": length, "minor_freq": minor_freq, "duration": res, "caller": caller, "repeat_type": repeat_type})

    return pl.DataFrame(records)


def test_no_interruption_str():
    minor_freq = 0.1
    motif_freqs = {"CAG": 1.0 - minor_freq, "CCG": minor_freq}
    return generate_motif_df(motif_freqs, minor_freq, "str_3")


def test_no_interruption_vntr(rng):
    minor_freq = 0.1
    motif_freqs = {random_dna(rng, 50): 1.0 - minor_freq, random_dna(rng, 50): minor_freq}
    return generate_motif_df(motif_freqs, minor_freq, "vntr_50")


def main() -> None:
    rng = np.random.default_rng(seed=1234)
    # Test cases:

    # Pure tandem repeats of different lengths (mostly a performance test / baseline), possibly split by different
    # lengths and/or motif sizes (STR vs VNTR)

    # Tandem repeats with different levels of alternate motif: 10%, 25%, 50%
    # TODO:
    res = pl.concat([test_no_interruption_str(), test_no_interruption_vntr(rng)])
    res.write_csv("./out/results.csv")

    # Tandem repeats with different levels of interruption: 5%, 10%, ... (although if TRviz doesn't support at all...)

    # Some kind of test of novel motifs...


if __name__ == "__main__":
    main()
