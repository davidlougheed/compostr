use parasail_rs::prelude::{Aligner, Alignment, Error, Matrix, Table};
use std::cmp;
use std::str::{self, Utf8Error};

use crate::motif::MotifSet;

/// Structure representing a computed motif decomposition, the result of a call to MotifSequenceDecomposer.decompose().
/// Contains the decomposition of the sequence into canonical motifs/sequence chunks + a CIGAR alignment for each
/// decomposed element (TODO) + the final score of the decomposition (i.e., interval-schedule weight) + TODO...
pub struct MotifSequenceDecomposition {
    pub decomposition: Vec<MotifAlignmentInterval>,
    pub score: i32, // Total weight achieved
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
    match_score: i32,
    mismatch_score: i32,
    gap_penalty: i32,
    aligner: Aligner,
    // interval computation parameters
    motif_alignment_score_cutoff: Option<i32>,
}

/// Representation of a motif alignment to a sequence.
/// Format: start (inclusive 0-based), end (inclusive 0-based), score, alignment table index
#[derive(Clone, Debug)]
pub struct MotifAlignmentInterval {
    start: usize,
    end: usize,
    score: i32,
    cigar: String,
    motif_idx: usize,
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
    if intervals[j - 1].score + m[p[j]] >= m[j - 1] {
        let temp: MotifAlignmentInterval = intervals[j - 1].clone();
        final_schedule.push(temp);
        backtrack_schedule(final_schedule, intervals, m, p, p[j])
    } else {
        backtrack_schedule(final_schedule, intervals, m, p, j - 1)
    }
}

/// Implementation of known weighted interval scheduling algorithm to do the motif decomposition
/// See https://en.wikipedia.org/wiki/Interval_scheduling#Weighted
/// Returns the best (or one of the best) schedules + its score.
fn schedule(intervals: &[MotifAlignmentInterval]) -> (Vec<MotifAlignmentInterval>, i32) {
    // sort intervals globally by earliest to latest end index
    let mut s_intervals: Vec<&MotifAlignmentInterval> = intervals.iter().collect();
    s_intervals.sort_by(|x, y| x.end.cmp(&y.end));

    // p[j] is the index of the latest interval that ends before interval j begins
    let mut p = vec![0; s_intervals.len() + 1];
    for i in 1..s_intervals.len() + 1 {
        let mut n = i - 1;
        while n > 0 {
            if s_intervals[n].end < s_intervals[i - 1].start {
                // interval coordinates are inclusive, so use '<'
                p[i] = n + 1;
                break;
            }
            n -= 1;
        }
    }

    // construct score table
    let mut m = vec![0; s_intervals.len() + 1];
    for i in 1..s_intervals.len() + 1 {
        let a = s_intervals[i - 1].score + m[p[i]];
        m[i] = cmp::max(a, m[i - 1]);
    }

    let mut final_schedule: Vec<MotifAlignmentInterval> = Vec::new();

    backtrack_schedule(&mut final_schedule, &s_intervals, &m, &p, s_intervals.len());

    final_schedule.reverse();

    (final_schedule, m[s_intervals.len()])
}

/// Encodes a CIGAR-ish alignment operation (insertion/deletion/match/mismatch).
/// With a grouping function, this can become a real CIGAR.
#[derive(Clone, Debug, PartialEq)]
enum AlignmentItem {
    Ins,
    Del,
    Match,
    Mismatch,
}

impl AlignmentItem {
    fn to_cigar_char(&self) -> char {
        match self {
            AlignmentItem::Ins => 'I',
            AlignmentItem::Del => 'D',
            AlignmentItem::Match => '=',
            AlignmentItem::Mismatch => 'X',
        }
    }
}

fn alignment_items_to_cigar(items: &[AlignmentItem]) -> String {
    let mut cigar = String::new();
    let mut ii = items.iter();
    if let Some(first) = ii.next() {
        let mut current_op = first;
        let mut current_count: usize = 1;
        for op in ii {
            if op == current_op {
                current_count += 1;
            } else {
                cigar = format!("{}{}{}", cigar, current_count, current_op.to_cigar_char());
                current_op = op;
                current_count = 1;
            }
        }
        cigar = format!("{}{}{}", cigar, current_count, current_op.to_cigar_char());
    }
    cigar
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
            .allow_ref_gaps(vec![String::from("prefix"), String::from("suffix")])
            .use_table()
            .striped()
            .build();

        MotifSequenceDecomposer {
            motif_set,
            match_score,
            mismatch_score,
            gap_penalty,
            aligner,
            motif_alignment_score_cutoff,
        }
    }

    /// Given a motif-sequence alignment table and an ending row/col for an alignment, trace back the alignment to return
    /// an interval tuple of: (starting sequence position, ending sequence position, score, CIGAR for the interval).
    /// This tuple will be used to build a MotifAlignmentInterval struct downstream.
    fn get_interval_from_score_matrix_start_pos(
        &self,
        seq: &[u8],
        motif: &[u8],
        tbl: &Table,
        mut row: usize,
        end_col: usize,
        cutoff: i32,
    ) -> Option<(usize, usize, i32, String)> {
        let mut col = end_col;

        let score = tbl.get(row, col);

        if let Some(s) = score
            && s >= cutoff
        {
            // keep alignment of motif to sequence from traceback as well:
            let mut motif_alignment: Vec<AlignmentItem> = Vec::new();

            while row > 0 {
                let mut options: Vec<(usize, usize, i32, AlignmentItem)> = Vec::new();
                if col > 0 {
                    // TODO: this doesn't support affine gap properly
                    if let Some(left) = tbl.get(row, col - 1) {
                        options.push((row, col - 1, left - self.gap_penalty, AlignmentItem::Ins));
                    }
                    if let Some(diag) = tbl.get(row - 1, col - 1) {
                        let (sc, ait) = if motif[row - 1] == seq[col - 1] {
                            (self.match_score, AlignmentItem::Match)
                        } else {
                            (self.mismatch_score, AlignmentItem::Mismatch)
                        };
                        options.push((row - 1, col - 1, diag + sc, ait));
                    }
                }
                // TODO: this doesn't support affine gap properly
                if let Some(up) = tbl.get(row - 1, col) {
                    options.push((row - 1, col, up - self.gap_penalty, AlignmentItem::Del));
                }
                let maxopt = options.iter().reduce(|acc, opt| if opt.2 > acc.2 { opt } else { acc });

                // maxopt shouldn't ever actually be None, otherwise something went wrong with score retrieval somehow.
                if let Some(mo) = maxopt {
                    row = mo.0;
                    col = mo.1;
                    motif_alignment.push(mo.3.clone());
                } else {
                    return None; // Something went wrong with score retrieval, this shouldn't happen
                }
            }

            // fixup cigar - we are always starting with a match or mismatch, because we get gaps at the front for free.
            //  TODO: validate this + maybe we want to be able to have bases at the front in the seq that become part of
            //   the motif?
            motif_alignment.push(if motif[row] == seq[col] { AlignmentItem::Match } else { AlignmentItem::Mismatch });
            motif_alignment.reverse();
            let cigar = alignment_items_to_cigar(&motif_alignment);

            return Some((col, end_col, s, cigar));
        }

        None
    }

    /// Given a set of motif alignments (with scoring tables) from Parasail plus an optional scoring cutoff, this function
    /// creates a vector of possible motif alignment intervals in the sequence. These will then be "scheduled" to produce
    /// the sequence motif decomposition.
    fn compute_intervals(
        &self,
        seq: &[u8],
        alignments: &[(&[u8], Alignment)],
        motif_alignment_score_cutoff: Option<i32>,
        seq_len: usize,
    ) -> Result<Vec<MotifAlignmentInterval>, Error> {
        let mut intervals: Vec<MotifAlignmentInterval> = Vec::with_capacity(alignments.len());
        let cutoff = motif_alignment_score_cutoff.unwrap_or(i32::MIN);
        for (ai, f) in alignments.iter().enumerate() {
            let motif_size = f.0.len() - 1;
            let tbl = f.1.get_score_table()?;
            eprintln!("scr {}", tbl);

            intervals.extend((0..seq_len).filter_map(|i| {
                self.get_interval_from_score_matrix_start_pos(seq, f.0, &tbl, motif_size, i, cutoff)
                    .map(|iv| MotifAlignmentInterval {
                        start: iv.0,
                        end: iv.1,
                        score: iv.2,
                        cigar: iv.3,
                        motif_idx: ai,
                    })
            }));
        }
        Ok(intervals)
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
        let intervals = self.compute_intervals(seq, &alignments, self.motif_alignment_score_cutoff, seq.len())?;
        eprintln!("{:?}", intervals);

        //  3: use weighted interval scheduling algorithm https://en.wikipedia.org/wiki/Interval_scheduling#Weighted
        //     to find best sequence of motifs, with any 'idle' time being non-motif DNA in between motifs.
        let (decomposition, score) = schedule(&intervals);

        Ok(MotifSequenceDecomposition { decomposition, score })
    }

    pub fn decomp_to_str(&self, decomp: MotifSequenceDecomposition) -> Result<Vec<&str>, Utf8Error> {
        let mut return_vec = Vec::new();
        for i in decomp.decomposition.iter() {
            let decomp_str: &str = str::from_utf8(self.motif_set.motifs[i.motif_idx].as_slice())?;
            return_vec.push(decomp_str);
        }
        Ok(return_vec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(b"CAGCAGCAGCAGCAGCAGCAGCAGCAG".to_vec(), vec!["CAG", "CAG", "CAG", "CAG", "CAG", "CAG", "CAG", "CAG", "CAG"])]
    #[case(b"CCGCCGCCGCCGCCGCCGCCGCCGCCG".to_vec(), vec!["CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG"])]
    #[case(b"CAGCAGCGGCAGCAAG".to_vec(), vec!["CAG", "CAG", "CGG", "CAG", "CAAG"])]
    #[case(b"CAGCAGCAAGTTCAGCCGCCGCCCG".to_vec(), vec!["CAG", "CAG", "CAAG", "T", "T", "CAG", "CCG", "CCG", "CCCG"])]
    fn test_decomposition(#[case] seq: Vec<u8>, #[case] expected_decomp: Vec<&str>) {
        let motif_set = MotifSet::new_from_strs(&vec!["CAG", "CCG"]);
        let decomposer = MotifSequenceDecomposer::new(motif_set, 5, -7, 4, Some(1));
        let res = decomposer.decompose(seq.as_slice()).unwrap();
        assert_eq!(decomposer.decomp_to_str(res).unwrap(), expected_decomp);
    }
}
