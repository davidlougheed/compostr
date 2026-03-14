use std::str::{self, Utf8Error};
use parasail_rs::prelude::{Aligner, Alignment, Error, Matrix};

use crate::motif::MotifSet;

/// TODO
pub struct MotifSequenceDecomposition {
    pub decomposition: Vec<Vec<u8>>,
    pub score: i32,  // Total weight achieved
}

impl MotifSequenceDecomposition {
    pub fn decomposition_strs<'m>(&'m self) -> Result<Vec<&'m str>, Utf8Error> {
        let mut res = Vec::with_capacity(self.decomposition.len());
        for m in self.decomposition.iter() {
            res.push(str::from_utf8(m)?);
        }
        Ok(res)
    }
}

/// 'Machine' that decomposes repetitive, TR-esque DNA sequences using a set of provided motifs.
pub struct MotifSequenceDecomposer {
    // scoring parameters:  TODO: more advanced scoring/matrix
    match_score: i32,
    mismatch_score: i32,
    gap_penalty: i32,
    // set of 'canonical' motifs for decomposition:
    pub motif_set: MotifSet,
}

/// Given a set of motif alignments (with scoring tables) from Parasail plus an optional scoring cutoff, this function
/// creates a vector of possible motif alignment intervals in the sequence. These will then be "scheduled" to produce
/// the sequence motif decomposition.
fn compute_intervals(
    alignments: &Vec<(&[u8], Alignment)>,
    motif_alignment_score_cutoff: Option<i32>,
    seq_len: usize,
) -> Result<Vec<(usize, usize, i32)>, Error> {
    let mut intervals: Vec<(usize, usize, i32)> = Vec::with_capacity(alignments.len());
    for f in alignments.iter() {
        let tbl = f.1.get_score_table()?;
        eprintln!("{}", tbl);

        for i in 0..seq_len {
            let s = tbl.get(f.0.len() - 1, i).expect("score table entry must exist");
            if let Some(cutoff) = motif_alignment_score_cutoff && s < cutoff { continue; }
            // TODO: do traceback
            intervals.push((0, i, s));  // TODO: not 0 but from traceback
        }
    }
    Ok(intervals)
}

/// TODO
fn schedule(alignments: &Vec<(&[u8], Alignment)>, intervals: &Vec<(usize, usize, i32)>) -> (Vec<Vec<u8>>, i32) {
    // re-use known weighted interval scheduling algorithm to do the motif decomposition
    // https://en.wikipedia.org/wiki/Interval_scheduling#Weighted

    // TODO

    // build up a decomposition of motifs, our "schedule"
    let decomposition = Vec::new();

    // TODO

    (decomposition, 0i32)  // TODO: real score
}

impl MotifSequenceDecomposer {
    pub fn new(motif_set: MotifSet, match_score: i32, mismatch_score: i32, gap_penalty: i32) -> Self {
        MotifSequenceDecomposer { motif_set, match_score, mismatch_score, gap_penalty }
    }

    /// TODO
    pub fn decompose<'m>(&'m self, seq: &[u8]) -> Result<MotifSequenceDecomposition, Error> {
        let motif_alignment_score_cutoff = Some(0);

        // rough algorithm outline, 3 parts:

        //  1: align all motifs (ends-free) to sequence to get alignment score matrix
        let matrix = Matrix::create(b"ACGT", self.match_score, self.mismatch_score).unwrap();
        let aligner = Aligner::new()
            .matrix(matrix)
            .gap_open(self.gap_penalty)
            .gap_extend(self.gap_penalty)
            .semi_global()
            .use_table()
            .striped()
            .build();

        let mut alignments: Vec<(&[u8], Alignment)> = Vec::with_capacity(self.motif_set.motifs.len());
        for m in self.motif_set.motifs.iter() {
            alignments.push((m, aligner.align(Some(m), seq)?));
        }

        //  2. determine intervals using some kind of heuristic so we don't have an absurd number?
        //     or just use last row(?) of the matrix as the score + figure out the interval... + do a little trimming
        let intervals = compute_intervals(&alignments, motif_alignment_score_cutoff, seq.len())?;

        //  3: use weighted interval scheduling algorithm https://en.wikipedia.org/wiki/Interval_scheduling#Weighted
        //     to find best sequence of motifs, with any 'idle' time being non-motif DNA in between motifs.
        let (decomposition, score) = schedule(&alignments, &intervals);

        Ok(MotifSequenceDecomposition { decomposition, score })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decomposition() {
        let motif_set = MotifSet::new_from_strs(&vec!["CAG", "CCG"]);
        let decomposer = MotifSequenceDecomposer::new(motif_set, 2, -7, 5);

        let res1 = decomposer.decompose(b"CAGCAGCAAGTTCAGCCGCCGCCCG").unwrap();
        assert_eq!(
            res1.decomposition_strs().unwrap(),
            vec!["CAG", "CAG", "CAAG", "T", "T", "CAG", "CCG", "CCG", "CCCG"]
        );
    }
}
