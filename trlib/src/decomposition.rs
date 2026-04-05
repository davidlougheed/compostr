use parasail_rs::prelude::{Aligner, Alignment, Error, Table};
use smallvec::SmallVec;
use std::cmp;
use std::sync::Arc;

use crate::motif::MotifSet;
use crate::scoring::{ScoringMatrix, ScoringMatrixError};

/// Structure representing a computed motif decomposition, the result of a call to MotifSequenceDecomposer.decompose().
/// Contains the decomposition of the sequence into canonical motifs/sequence chunks + a CIGAR alignment for each
/// decomposed element (TODO) + the final score of the decomposition (i.e., interval-schedule weight) + TODO...
pub struct MotifSequenceDecomposition {
    pub motif_set: Arc<MotifSet>,
    pub decomposition: Vec<DecompositionItem>,
    pub score: i32,    // Total weight achieved
    pub copies: usize, // Total number of copies of any motif
}

pub enum DecompositionItem {
    Alignment(MotifAlignmentInterval),
    Gap(usize),
}

impl MotifSequenceDecomposition {
    /// Returns a vector of CIGAR items representing the alignment of a sequence just composed of motifs from the
    /// motif set against the original sequence.
    pub fn cigar(&self) -> Vec<CigarItem> {
        let mut cigars = Vec::<CigarItem>::new();
        for d in self.decomposition.iter() {
            match d {
                DecompositionItem::Alignment(a) => {
                    cigars.extend_from_slice(&a.cigar);
                }
                DecompositionItem::Gap(g) => {
                    // Infill any gaps where we didn't decompose to any motif at all
                    cigars.push(CigarItem::Ins(*g));
                }
            }
        }
        cigars
    }

    /// Returns a string representation of .cigar()
    pub fn cigar_string(&self) -> String {
        let cigar_strings: Vec<String> = self.cigar().into_iter().map(|c| c.to_alignment_string()).collect();
        cigar_strings.join("")
    }

    /// Given the original sequence, returns the decomposition of it into sub-slices.
    pub fn sequence_items<'a>(&self, seq: &'a [u8], with_unmapped: bool) -> Vec<&'a [u8]> {
        let mut res: Vec<&'a [u8]> = Vec::new();
        let mut last_end = 0;
        for d in self.decomposition.iter() {
            match d {
                DecompositionItem::Alignment(a) => {
                    res.push(&seq[a.start..a.end + 1]);
                    last_end = a.end;
                }
                DecompositionItem::Gap(g) => {
                    if with_unmapped {
                        res.push(&seq[last_end + 1..last_end + g + 1]);
                    }
                }
            };
        }
        res
    }

    /// Returns a human-readable alignment string representation given the originally decomposed string, which will
    /// look something like:
    ///
    /// CAGCAGCAGCAGCA-G
    /// |||||||X|||||| |
    /// CAGCAGCGGCAGCAAG
    ///
    /// Where the top string is generated from the motif set + the found decomposition, and the bottom string is the
    /// original sequence.
    pub fn alignment_string(&self, seq: &[u8]) -> String {
        let mut query_strings = Vec::<String>::new();
        let mut align_strings = Vec::<String>::new();
        let mut seq_strings = Vec::<String>::new();

        let mut last_end = 0;
        for d in self.decomposition.iter() {
            match d {
                DecompositionItem::Alignment(a) => {
                    let m = &self.motif_set.motifs[a.motif_idx];
                    let mut qi = 0;
                    let mut ri = a.start;
                    for item in a.cigar.iter() {
                        align_strings.push(item.to_alignment_string());
                        match item {
                            CigarItem::Del(c) => {
                                query_strings.push(String::from_utf8_lossy(&m[qi..qi + c]).to_string());
                                seq_strings.push("-".repeat(*c));
                                qi += c;
                            }
                            CigarItem::Ins(c) => {
                                query_strings.push("-".repeat(*c));
                                seq_strings.push(String::from_utf8_lossy(&seq[ri..ri + c]).to_string());
                                ri += c;
                            }
                            CigarItem::Match(c) | CigarItem::Mismatch(c) => {
                                query_strings.push(String::from_utf8_lossy(&m[qi..qi + c]).to_string());
                                seq_strings.push(String::from_utf8_lossy(&seq[ri..ri + c]).to_string());
                                qi += c;
                                ri += c;
                            }
                        }
                    }
                    last_end = a.end;
                }
                DecompositionItem::Gap(g) => {
                    // Infill any gaps where we didn't decompose to any motif at all
                    let undecomp_seq = String::from_utf8_lossy(&seq[last_end + 1..last_end + g + 1]).to_string();
                    query_strings.push(" ".repeat(undecomp_seq.len())); // different type of gap, use space
                    align_strings.push(" ".repeat(undecomp_seq.len()));
                    seq_strings.push(undecomp_seq);
                }
            }
        }

        format!(
            "{}\n{}\n{}",
            query_strings.join(""),
            align_strings.join(""),
            seq_strings.join("")
        )
    }
}

/// 'Machine' that decomposes repetitive, TR-esque DNA sequences using a set of provided motifs.
pub struct MotifSequenceDecomposer {
    // set of 'canonical' motifs for decomposition:
    pub motif_set: Arc<MotifSet>,
    // alignment:  TODO: more advanced scoring/matrix
    match_score: i32,
    mismatch_score: i32,
    gap_penalty: i32,
    aligner: Aligner,
    // interval computation parameters
    motif_alignment_score_cutoff: i32,
}

/// Representation of a motif alignment to a sequence.
#[derive(Clone, Debug)]
pub struct MotifAlignmentInterval {
    start: usize, // inclusive, 0-based
    end: usize,   // inclusive, 0-based
    pub score: i32,
    pub cigar: SmallVec<[CigarItem; 4]>,
    motif_idx: usize, // Index of the motif in the motif set
}

fn backtrack_schedule(
    final_schedule: &mut Vec<MotifAlignmentInterval>,
    mut intervals: Vec<MotifAlignmentInterval>,
    m: &[i32],
    p: &[usize],
    mut j: usize,
) {
    while j > 0 {
        if intervals[j - 1].score + m[p[j]] >= m[j - 1] {
            // avoid clone by moving intervals[j - 1]; we're done with everything from [j - 1] forward anyway
            // since we only ever decrement j (via p or via subtracting 1).
            let next_item = {
                let mut d = intervals.drain((j - 1)..);
                d.next().expect("backtract_schedule: next item must exist")
            };
            final_schedule.push(next_item);
            j = p[j];
        } else {
            j -= 1;
        }
    }
}

/// Implementation of known weighted interval scheduling algorithm to do the motif decomposition
/// See https://en.wikipedia.org/wiki/Interval_scheduling#Weighted
/// Requires an input vector already sorted by end position.
/// Returns the best (or one of the best) schedules + its score.
fn schedule(intervals: Vec<MotifAlignmentInterval>) -> (Vec<MotifAlignmentInterval>, i32) {
    // p[j] is the index of the latest interval that ends before interval j begins
    let mut p = vec![0; intervals.len() + 1];
    for i in 1..intervals.len() + 1 {
        let mut n = i - 1;
        while n > 0 {
            if intervals[n].end < intervals[i - 1].start {
                // interval coordinates are inclusive, so use '<'
                p[i] = n + 1;
                break;
            }
            n -= 1;
        }
    }

    // construct score table
    let mut m = vec![0; intervals.len() + 1];
    for i in 1..intervals.len() + 1 {
        let a = intervals[i - 1].score + m[p[i]];
        m[i] = cmp::max(a, m[i - 1]);
    }

    let mut final_schedule: Vec<MotifAlignmentInterval> = Vec::new();
    let n_intervals = intervals.len();

    backtrack_schedule(&mut final_schedule, intervals, &m, &p, n_intervals);

    final_schedule.reverse();

    (final_schedule, m[n_intervals])
}

/// Given a schedule of MotifAlignmentInterval, wrap with DecompositionItem::Alignment and fill in gaps with
/// DecompositionItem::Gap so that we have a contiguous record of alignment/gap for the whole sequence.
fn interval_schedule_to_decomposition(schedule: Vec<MotifAlignmentInterval>, seq_len: usize) -> Vec<DecompositionItem> {
    let mut last_end = 0;
    let mut decomposition = Vec::with_capacity(schedule.len()); // Size >= schedule.len()
    for d in schedule.into_iter() {
        if d.start > last_end + 1 {
            // Infill any gaps where we didn't decompose to any motif at all
            decomposition.push(DecompositionItem::Gap(d.start - last_end - 1));
        }
        last_end = d.end;
        decomposition.push(DecompositionItem::Alignment(d));
    }
    if last_end < seq_len - 1 {
        decomposition.push(DecompositionItem::Gap(seq_len - last_end - 1));
    }
    decomposition
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
    fn to_cigar_item(&self, count: usize) -> CigarItem {
        match self {
            AlignmentItem::Ins => CigarItem::Ins(count),
            AlignmentItem::Del => CigarItem::Del(count),
            AlignmentItem::Match => CigarItem::Match(count),
            AlignmentItem::Mismatch => CigarItem::Mismatch(count),
        }
    }
}

#[derive(Clone, Debug)]
pub enum CigarItem {
    Ins(usize),
    Del(usize),
    Match(usize),
    Mismatch(usize),
}

impl CigarItem {
    pub fn to_cigar_string(&self) -> String {
        match self {
            Self::Ins(c) => format!("{c}I"),
            Self::Del(c) => format!("{c}D"),
            Self::Match(c) => format!("{c}="),
            Self::Mismatch(c) => format!("{c}X"),
        }
    }

    fn to_alignment_string(&self) -> String {
        match self {
            Self::Ins(c) | Self::Del(c) => " ".repeat(*c),
            Self::Match(c) => "|".repeat(*c),
            Self::Mismatch(c) => "X".repeat(*c),
        }
    }
}

impl MotifSequenceDecomposer {
    pub fn new(
        motif_set: MotifSet,
        match_score: i32,
        mismatch_score: i32,
        gap_penalty: i32,
        motif_alignment_score_cutoff: Option<i32>,
    ) -> Result<Self, ScoringMatrixError> {
        let matrix = ScoringMatrix::new_iupac_dna(match_score, mismatch_score)?;
        let aligner = Aligner::new()
            .matrix(matrix.matrix)
            .gap_open(gap_penalty)
            .gap_extend(gap_penalty)
            .semi_global()
            .allow_ref_gaps(vec![String::from("prefix"), String::from("suffix")])
            .use_table()
            .striped()
            .build();

        Ok(MotifSequenceDecomposer {
            motif_set: Arc::new(motif_set),
            match_score,
            mismatch_score,
            gap_penalty,
            aligner,
            motif_alignment_score_cutoff: motif_alignment_score_cutoff.unwrap_or(i32::MIN),
        })
    }

    /// Maximum-scoring option function for semi-global alignment traceback
    fn sg_traceback_max_opt(
        &self,
        motif: &[u8],
        seq: &[u8],
        row: usize,
        col: usize,
        tbl_slice: &[i32],
        tbl_cols: usize,
    ) -> Option<(usize, usize, i32, AlignmentItem)> {
        // Instead of keeping options in vec, save a lot of time by just enumerating every possible comparison
        let (ins_opt, match_opt) = if col > 0 {
            // TODO: this doesn't support affine gap properly
            let left = tbl_slice[row * tbl_cols + (col - 1)];
            let diag = tbl_slice[(row - 1) * tbl_cols + (col - 1)];
            let (sc, ait) = if motif[row - 1] == seq[col - 1] {
                (self.match_score, AlignmentItem::Match)
            } else {
                (self.mismatch_score, AlignmentItem::Mismatch)
            };
            (
                // insertion
                Some((row, col - 1, left - self.gap_penalty, AlignmentItem::Ins)),
                // match or mismatch
                Some((row - 1, col - 1, diag + sc, ait)),
            )
        } else {
            (None, None)
        };
        // TODO: this doesn't support affine gap properly
        let up = tbl_slice[(row - 1) * tbl_cols + col];
        let del_opt = (row - 1, col, up - self.gap_penalty, AlignmentItem::Del);

        let maxopt = match (ins_opt, match_opt) {
            (Some(i), Some(m)) => {
                Some(if i.2 > m.2 && i.2 >= del_opt.2 {
                    i
                } else if m.2 >= i.2 && m.2 >= del_opt.2 {
                    m
                } else {
                    // if d.2 > m.2 && d.2 > i.2
                    del_opt
                })
            }
            (Some(i), None) => Some(if i.2 >= del_opt.2 { i } else { del_opt }),
            (None, Some(m)) => Some(if m.2 >= del_opt.2 { m } else { del_opt }),
            // base case (should never happen)
            (None, None) => None,
        };

        maxopt
    }

    /// Given a motif-sequence alignment table and an ending row/col for an alignment, trace back the alignment to return
    /// an interval tuple of: (starting sequence position, ending sequence position, score, CIGAR for the interval).
    /// This tuple will be used to build a MotifAlignmentInterval struct downstream.
    fn get_interval_from_score_matrix_start_pos(
        &self,
        seq: &[u8],
        motif: &[u8],
        motif_idx: usize,
        tbl: &Table,
        mut row: usize,
        end_col: usize,
    ) -> Option<MotifAlignmentInterval> {
        let mut col = end_col;
        let tbl_slice = tbl.as_slice();
        let tbl_cols = tbl.cols();

        let score = tbl_slice[row * tbl_cols + col];

        if score >= self.motif_alignment_score_cutoff {
            // keep alignment of motif to sequence from traceback as well:
            let mut current_op: AlignmentItem = AlignmentItem::Match; // Dummy value to be replaced
            let mut current_op_count: usize = 0;

            let mut cigar: SmallVec<[CigarItem; 4]> = SmallVec::new();

            while row > 0 {
                let maxopt = self.sg_traceback_max_opt(motif, seq, row, col, tbl_slice, tbl_cols);

                // maxopt shouldn't ever actually be None, otherwise something went wrong with score retrieval somehow.
                if let Some(mo) = maxopt {
                    row = mo.0;
                    col = mo.1;
                    if mo.3 != current_op {
                        if current_op_count > 0 {
                            cigar.push(current_op.to_cigar_item(current_op_count));
                        }
                        current_op = mo.3;
                        current_op_count = 1;
                    } else {
                        current_op_count += 1;
                    }
                } else {
                    return None; // Something went wrong with score retrieval, this shouldn't happen
                }
            }

            // fixup cigar - we are always starting with a match or mismatch, because we get gaps at the front for free.
            //  TODO: validate this + maybe we want to be able to have bases at the front in the seq that become part of
            //   the motif?
            let last = if motif[row] == seq[col] {
                AlignmentItem::Match
            } else {
                AlignmentItem::Mismatch
            };

            if last != current_op {
                if current_op_count > 0 {
                    cigar.push(current_op.to_cigar_item(current_op_count));
                }
                current_op = last;
                current_op_count = 1;
            } else {
                current_op_count += 1;
            }

            cigar.push(current_op.to_cigar_item(current_op_count));
            cigar.reverse();

            return Some(MotifAlignmentInterval {
                start: col,
                end: end_col,
                score,
                cigar,
                motif_idx,
            });
        }

        None
    }

    /// Given a set of motif alignments (with scoring tables) from Parasail plus an optional scoring cutoff, this function
    /// creates a vector of possible motif alignment intervals in the sequence. These will then be "scheduled" to produce
    /// the sequence motif decomposition. The returned vector is naturally already sorted by end position.
    fn compute_intervals(
        &self,
        seq: &[u8],
        alignments: Vec<(&[u8], Alignment)>,
    ) -> Result<Vec<MotifAlignmentInterval>, Error> {
        let mut intervals: Vec<MotifAlignmentInterval> = Vec::with_capacity(alignments.len() * seq.len());

        let mut tables: Vec<Table> = Vec::with_capacity(alignments.len());
        for (_, a) in alignments.iter() {
            tables.push(a.get_score_table()?);
        }

        for i in 0..seq.len() {
            let mut best_intervals_for_pos = SmallVec::<[MotifAlignmentInterval; 2]>::new();
            let mut best_score: i32 = i32::MIN;
            let mut best_score_len: usize = usize::MAX;

            for ai in 0..alignments.len() {
                let motif = alignments[ai].0;
                let tbl = &tables[ai];
                let iv = self.get_interval_from_score_matrix_start_pos(seq, motif, ai, tbl, motif.len() - 1, i);
                if let Some(interval) = iv {
                    if (interval.score > best_score && interval.end - interval.start <= best_score_len)
                        || (interval.score == best_score && interval.end - interval.start < best_score_len)
                    {
                        // cases where the interval is strictly better than the best one found so far:
                        //  - it is above the best score and at most the same length (will replace in schedule for a better result)
                        //  - it is the same score and shorter (will replace in schedule for more unallocated time)
                        best_score = interval.score;
                        best_score_len = interval.end - interval.start;
                        best_intervals_for_pos.clear();
                        best_intervals_for_pos.push(interval);
                    } else if interval.score == best_score && interval.end - interval.start == best_score_len {
                        // case where the interval is the same (for future use outputting multiple schedules):
                        best_intervals_for_pos.push(interval);
                    }
                }
            }

            intervals.extend(best_intervals_for_pos);
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

        //  2. determine intervals, cutting off low-scoring possibilities; returns a vector already sorted by end pos.
        let intervals = self.compute_intervals(seq, alignments)?;

        //  3: use weighted interval scheduling algorithm https://en.wikipedia.org/wiki/Interval_scheduling#Weighted
        //     to find best sequence of motifs, with any 'idle' time being non-motif DNA in between motifs.
        let (schedule, score) = schedule(intervals);
        let copies = schedule.len(); // Copy number of tandem repeat

        //  4: build final decomposition items by interspersing gaps as needed
        let decomposition = interval_schedule_to_decomposition(schedule, seq.len());

        Ok(MotifSequenceDecomposition {
            motif_set: self.motif_set.clone(),
            decomposition,
            score,
            copies,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(b"CAG".to_vec(), vec!["CAG"], "CAG\n|||\nCAG", 1)]
    #[case(b"CAAG".to_vec(), vec!["CAAG"], "CA-G\n|| |\nCAAG", 1)]
    #[case(
        b"CAGCAGCAGCAGCAGCAGCAGCAGCAG".to_vec(),
        vec!["CAG", "CAG", "CAG", "CAG", "CAG", "CAG", "CAG", "CAG", "CAG"],
        "CAGCAGCAGCAGCAGCAGCAGCAGCAG\n\
         |||||||||||||||||||||||||||\n\
         CAGCAGCAGCAGCAGCAGCAGCAGCAG",
        9,
    )]
    #[case(
        b"CCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCAG".to_vec(),
        vec!["CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG",
             "CAG"],
        "CCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCAG\n\
         ||||||||||||||||||||||||||||||||||||||||||||||||\n\
         CCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCAG",
        16,
    )]
    #[case(
        b"CAGCAGCGGCAGCAAG".to_vec(),
        vec!["CAG", "CAG", "CGG", "CAG", "CAAG"],
        "CAGCAGCCGCAGCA-G\n\
         |||||||X|||||| |\n\
         CAGCAGCGGCAGCAAG",
        5,
    )]
    #[case(
        b"CAGCAGCAAGTTCAGCCGCCGCCCG".to_vec(),
        vec!["CAG", "CAG", "CAAG", "TT", "CAG", "CCG", "CCG", "CCCG"],
        "CAGCAGCA-G  CAGCCGCCGCC-G\n\
         |||||||| |  ||||||||||| |\n\
         CAGCAGCAAGTTCAGCCGCCGCCCG",
        7,
    )]
    fn test_decomposition(
        #[case] seq: Vec<u8>,
        #[case] expected_decomp: Vec<&str>,
        #[case] expected_align_str: &str,
        #[case] expected_copies: usize,
    ) {
        let motif_set = MotifSet::new_from_strs(&vec!["CAG", "CCG"]);
        let decomposer = MotifSequenceDecomposer::new(motif_set, 5, -7, 4, Some(1)).unwrap();
        let res = decomposer.decompose(seq.as_slice()).unwrap();
        let expected_decomp_u8: Vec<&[u8]> = expected_decomp.iter().map(|c| c.as_bytes()).collect();
        assert_eq!(res.sequence_items(&seq, true), expected_decomp_u8);
        assert_eq!(res.alignment_string(&seq), expected_align_str);
        assert_eq!(res.copies, expected_copies);
    }

    // TODO: check this - the X are wrong in the actual result
    const L1: &str =
        "    GATGATGGGAGTGTGCGCAGTGTAAGGATGATGGGAGTGTGCGCAGTGT-AAGGATGATGGGAGTGTGCGCAGTGT-AAG                 \n";
    const L2: &str =
        "    |||||||||||||||||||||||||||||||||||||||||X|||X||| | ||||||||||||||||||X||||| | |                 \n";
    const L3: &str =
        "TGAGGATGATGGGAGTGTGCGCAGTGTAAGGATGATGGGAGTGTGTGCAATGTGA-GGATGATGGGAGTGTGCACAGTGTGA-GGACGATGGGAGTGTGCG";

    #[rstest]
    #[case(
        b"GTGAGGATGATGGGAGTGTGCGCAGTGTAAGGATGATGGGAGTGTGTGCAATGTGAGGATGATGGGAGTGTGCACAGTGTGAGGACGATGGGAGTGTGCG".to_vec(),
        &format!("{L1}{L2}{L3}"),
        3,
    )]
    fn test_decomposition_2(#[case] seq: Vec<u8>, #[case] expected_align_str: &str, #[case] expected_copies: usize) {
        let motif_set = MotifSet::new_from_strs(&vec!["GATGATGGGAGTGTGCGCAGTGTAAG"]);
        let decomposer = MotifSequenceDecomposer::new(motif_set, 5, -7, 4, Some(-2)).unwrap();
        let res = decomposer.decompose(seq.as_slice()).unwrap();
        assert_eq!(res.alignment_string(&seq), expected_align_str);
        assert_eq!(res.copies, expected_copies);
    }
}
