use std::str::{self, Utf8Error};
use parasail_rs::prelude::{Aligner, Alignment, Error, Matrix, Table};

use crate::motif::MotifSet;

/// Structure representing a computed motif decomposition, the result of a call to MotifSequenceDecomposer.decompose().
/// Contains the decomposition of the sequence into canonical motifs/sequence chunks + a CIGAR alignment for each
/// decomposed element (TODO) + the final score of the decomposition (i.e., interval-schedule weight) + TODO...
pub struct MotifSequenceDecomposition {
    pub decomposition: Vec<Vec<u8>>,  // Vector of string bytevectors right now, but this will probably change
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
    // set of 'canonical' motifs for decomposition:
    pub motif_set: MotifSet,
    // alignment:  TODO: more advanced scoring/matrix
    aligner: Aligner,
    // interval computation parameters
    motif_alignment_score_cutoff: Option<i32>,
}

fn get_interval_from_score_matrix_start_pos(
    tbl: &Table, start_row: usize, start_col: usize, cutoff: i32
) -> Option<(usize, usize, i32)> {
    let mut row = start_row;
    let mut col = start_col;

    let score = tbl.get(row, col);

    if let Some(s) = score && s >= cutoff {
        // TODO: keep alignment of motif to sequence from traceback as well.
        while row > 0 {
            let mut options: Vec<(usize, usize, i32)> = Vec::new();
            if col > 0 {
                if let Some(left) = tbl.get(row, col - 1) { options.push((row, col - 1, left)); }
                if let Some(diag) = tbl.get(row - 1, col - 1) { options.push((row - 1, col - 1, diag)); }
            }
            if let Some(up) = tbl.get(row - 1, col) { options.push((row - 1, col, up)); }
            let maxopt = options.iter().reduce(|acc, opt| if opt.2 > acc.2 { opt } else { acc } );

            // maxopt shouldn't ever actually be None, otherwise something went wrong with score retrieval somehow.
            if let Some(&mo) = maxopt {
                row = mo.0;
                col = mo.1;
            } else {
                return None  // Soemthing went wrong with score retrieval, this shouldn't happen
            }
        }

        return Some((col, start_col, s));
    }

    None
}

#[derive(Debug)]
struct MotifAlignmentInterval(usize, usize, i32, usize);

/// Given a set of motif alignments (with scoring tables) from Parasail plus an optional scoring cutoff, this function
/// creates a vector of possible motif alignment intervals in the sequence. These will then be "scheduled" to produce
/// the sequence motif decomposition.
fn compute_intervals(
    alignments: &Vec<(&[u8], Alignment)>,
    motif_alignment_score_cutoff: Option<i32>,
    seq_len: usize,
) -> Result<Vec<MotifAlignmentInterval>, Error> {
    let mut intervals: Vec<MotifAlignmentInterval> = Vec::with_capacity(alignments.len());
    let cutoff = motif_alignment_score_cutoff.unwrap_or(i32::MIN);
    for (ai, f) in alignments.iter().enumerate() {
        let motif_size = f.0.len() - 1;
        let tbl = f.1.get_score_table()?;
        eprintln!("{}", tbl);

        for i in 0..seq_len {
            if let Some(iv) = get_interval_from_score_matrix_start_pos(&tbl, motif_size, i, cutoff) {
                intervals.push(MotifAlignmentInterval(iv.0, iv.1, iv.2, ai));
            }
        }
    }
    Ok(intervals)
}

/// Implementation of known weighted interval scheduling algorithm to do the motif decomposition
/// See https://en.wikipedia.org/wiki/Interval_scheduling#Weighted
fn schedule(
    seq: &[u8],
    alignments: &Vec<(&[u8], Alignment)>,
    intervals: &Vec<MotifAlignmentInterval>, // Vector of tuples (start, end, score, alignment index)
) -> (Vec<Vec<u8>>, i32) {
    // TODO

    // build up a decomposition of motifs, our "schedule"
    let decomposition = Vec::new();

    // TODO

    (decomposition, 0i32)  // TODO: real score
}

impl MotifSequenceDecomposer {
    pub fn new(
        motif_set: MotifSet,
        match_score: i32,
        mismatch_score: i32,
        gap_penalty: i32,
        motif_alignment_score_cutoff: Option<i32>,
    ) -> Self {
        let matrix = Matrix::create(b"ACGT", match_score, mismatch_score).unwrap();
        let aligner = Aligner::new()
            .matrix(matrix)
            .gap_open(gap_penalty)
            .gap_extend(gap_penalty)
            .semi_global()
            .use_table()
            .striped()
            .build();

        MotifSequenceDecomposer { motif_set, aligner, motif_alignment_score_cutoff }
    }

    /// TODO
    pub fn decompose(&self, seq: &[u8]) -> Result<MotifSequenceDecomposition, Error> {
        // rough algorithm outline, 3 parts:

        //  1: align all motifs (ends-free) to sequence to get alignment score matrix
        let mut alignments: Vec<(&[u8], Alignment)> = Vec::with_capacity(self.motif_set.motifs.len());
        for m in self.motif_set.motifs.iter() {
            alignments.push((m, self.aligner.align(Some(m), seq)?));
        }

        //  2. determine intervals using some kind of heuristic so we don't have an absurd number?
        //     or just use last row(?) of the matrix as the score + figure out the interval... + do a little trimming
        let intervals = compute_intervals(&alignments, self.motif_alignment_score_cutoff, seq.len())?;
        eprintln!("{:?}", intervals);

        //  3: use weighted interval scheduling algorithm https://en.wikipedia.org/wiki/Interval_scheduling#Weighted
        //     to find best sequence of motifs, with any 'idle' time being non-motif DNA in between motifs.
        let (decomposition, score) = schedule(seq, &alignments, &intervals);

        Ok(MotifSequenceDecomposition { decomposition, score })
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use super::*;

    #[rstest]
    #[case(b"CAGCAGCAGCAGCAGCAGCAGCAGCAG".to_vec(), vec!["CAG", "CAG", "CAG", "CAG", "CAG", "CAG", "CAG", "CAG", "CAG"])]
    #[case(b"CAGCAGCGGCAGCAAG".to_vec(), vec!["CAG", "CAG", "CGG", "CAG", "CAAG"])]
    #[case(b"CAGCAGCAAGTTCAGCCGCCGCCCG".to_vec(), vec!["CAG", "CAG", "CAAG", "T", "T", "CAG", "CCG", "CCG", "CCCG"])]
    fn test_decomposition(#[case] seq: Vec<u8>, #[case] expected_decomp: Vec<&str>) {
        let motif_set = MotifSet::new_from_strs(&vec!["CAG", "CCG"]);
        let decomposer = MotifSequenceDecomposer::new(motif_set, 5, -7, 4, Some(1));
        let res = decomposer.decompose(seq.as_slice()).unwrap();
        assert_eq!(res.decomposition_strs().unwrap(), expected_decomp);
    }
}
