use std::str::{self, Utf8Error};
use parasail_rs::prelude::{Aligner, Alignment, Error, Matrix, Table};

use crate::motif::MotifSet;

/// Structure representing a computed motif decomposition, the result of a call to MotifSequenceDecomposer.decompose().
/// Contains the decomposition of the sequence into canonical motifs/sequence chunks + a CIGAR alignment for each
/// decomposed element (TODO) + the final score of the decomposition (i.e., interval-schedule weight) + TODO...
pub struct MotifSequenceDecomposition {
    pub decomposition: Vec<MotifAlignmentInterval>,
    pub score: i32,  // Total weight achieved
}

impl MotifSequenceDecomposition {
    pub fn decomposition_strs(&self) -> Result<Vec<&str>, Utf8Error> {
        let mut res = Vec::with_capacity(self.decomposition.len());
        for m in self.decomposition.iter() {
            //res.push(str::from_utf8(m.3)?);
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

/// Given a motif-sequence alignment table and an ending row/col for an alignment, trace back the alignment to return
/// an interval tuple of: (starting sequence position, ending sequence position, score)
fn get_interval_from_score_matrix_start_pos(
    tbl: &Table, mut row: usize, end_col: usize, cutoff: i32
) -> Option<(usize, usize, i32)> {
    let mut col = end_col;

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

        return Some((col, end_col, s));
    }

    None
}

/// Representation of a motif alignment to a sequence.
/// Format: start (inclusive 0-based), end (inclusive 0-based), score, alignment table index
#[derive(Debug)]
#[derive(Clone)]
struct MotifAlignmentInterval(usize, usize, i32, usize);

/// Given a set of motif alignments (with scoring tables) from Parasail plus an optional scoring cutoff, this function
/// creates a vector of possible motif alignment intervals in the sequence. These will then be "scheduled" to produce
/// the sequence motif decomposition.
fn compute_intervals(
    alignments: &[(&[u8], Alignment)],
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

fn backtrack_schedule(
    final_schedule: &mut Vec<MotifAlignmentInterval>,
    intervals: &Vec<&MotifAlignmentInterval>,
    m: &Vec<i32>,
    p: &Vec<usize>,
    j: usize,
) {
    if j == 0 {
        return;
    }
    if intervals[j-1].2 + m[p[j]] >= m[j - 1] {
        let temp: MotifAlignmentInterval = intervals[j-1].clone();
        final_schedule.push(temp);
        backtrack_schedule(final_schedule, intervals, m, p, p[j])
    } else {
        backtrack_schedule(final_schedule, intervals, m, p, j-1)
    }
}

/// Implementation of known weighted interval scheduling algorithm to do the motif decomposition
/// See https://en.wikipedia.org/wiki/Interval_scheduling#Weighted
fn schedule(
    seq: &[u8],
    alignments: &[(&[u8], Alignment)],
    intervals: &[MotifAlignmentInterval], // Vector of tuples (start, end, score, alignment index)
) -> (Vec<MotifAlignmentInterval>, i32) {

    //sort intervals globally by earliest to latest end index
    let mut s_intervals: Vec<&MotifAlignmentInterval> = Vec::new();
    for i in intervals.iter() {
        s_intervals.push(i);
    }
    s_intervals.sort_by(|x, y| x.1.cmp(&y.1));

    // p[j] is the index of the latest interval that ends before interval j begins
    let mut p = vec![0; s_intervals.len()+1];
    for i in 1..s_intervals.len() + 1 {
        let mut n = i - 1;
        while n > 0 {
            if s_intervals[n].1 < s_intervals[i-1].0 {
                p[i] = n + 1;
                break;
            }
            n -= 1;
        }
    }

    //construct score table
    let mut m = vec![0; s_intervals.len()+1];
    for i in 1..s_intervals.len() + 1 {
        let a = s_intervals[i-1].2 + m[p[i]];
        if a > m[i-1] {
            m[i] = a;
        } else {
            m[i] = m[i-1];
        }
    }

    let mut final_schedule: Vec<MotifAlignmentInterval> = Vec::new();

    backtrack_schedule(&mut final_schedule, &s_intervals, &m, &p, s_intervals.len());

    final_schedule.reverse();

    (final_schedule, m[s_intervals.len()])
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

    /// Main functionality for MotifSequenceDecomposer - given a sequence, decomposes it into motifs using the motif
    /// set and alignment parameters that were specified for the decomposer instance.
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

    pub fn decomp_to_str(&self, decomp: MotifSequenceDecomposition) -> Result<Vec<&str>, Utf8Error> {
        let mut return_vec = Vec::new();
        for i in decomp.decomposition.iter() {
            let decomp_str: &str = str::from_utf8(self.motif_set.motifs[i.3].as_slice())?;
            return_vec.push(decomp_str);
        }
        Ok(return_vec)
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
        assert_eq!(decomposer.decomp_to_str(res).unwrap(), expected_decomp);
    }
}
