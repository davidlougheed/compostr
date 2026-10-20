from __future__ import annotations

import numpy as np

BASES = "ACGT"
BASES_NP = np.array(list(BASES))


def random_dna(rng: np.random.Generator, length: int) -> str:
    return "".join(rng.choice(BASES_NP, size=length))


class RepeatUnit(str):
    """Representation of a repeat unit in a generated sequence"""


class Interruption(str):
    """Representation of an interruption in a generated sequence"""


class GeneratedSequence:
    def __init__(self, seq: list[RepeatUnit | Interruption]):
        self.seq = seq

    def append(self, item: RepeatUnit | Interruption):
        self.seq.append(item)

    def __str__(self):
        return "".join(self.seq)

    def __bytes__(self):
        return str(self).encode("ascii")


class Generator:
    def __init__(self, motif_probs: dict[str, float], p_interruption: float, interruption_params: tuple[int, float]):
        if sum(motif_probs.values()) != 1.0:
            raise ValueError("Motif probabilities must sum to 1.0")
        self.motifs = tuple(map(RepeatUnit, motif_probs.keys()))
        self.motif_probs = tuple(motif_probs.values())
        self.p_interruption = p_interruption
        self.neg_binom_r = interruption_params[0]
        self.neg_binom_p = interruption_params[1]

        self.rng = np.random.default_rng(seed=1234)

    def generate_sequence(self, n_motifs: int) -> GeneratedSequence:
        res = GeneratedSequence([])

        while n_motifs > 0:
            ru = self.rng.choice(self.motifs, p=self.motif_probs)
            res.append(ru)
            n_motifs -= 1
            if self.p_interruption and self.rng.random() < self.p_interruption:
                interruption_size = self.rng.negative_binomial(self.neg_binom_r, self.neg_binom_p)
                interruption = "".join(self.rng.choice(BASES_NP, size=interruption_size))
                res.append(Interruption(interruption))

        return res
