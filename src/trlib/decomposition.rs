use parasail_rs::prelude::{Aligner, Alignment, Error, Matrix, Table};
use std::cmp;
use std::str::{self, Utf8Error};
use std::sync::Arc;

use crate::motif::MotifSet;

/// Structure representing a computed motif decomposition, the result of a call to MotifSequenceDecomposer.decompose().
/// Contains the decomposition of the sequence into canonical motifs/sequence chunks + a CIGAR alignment for each
/// decomposed element (TODO) + the final score of the decomposition (i.e., interval-schedule weight) + TODO...
pub struct MotifSequenceDecomposition {
    pub motif_set: Arc<MotifSet>,
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

    pub fn cigar(&self) -> Vec<CigarItem> {
        self.decomposition.iter().flat_map(|d| d.cigar.clone()).collect()
    }

    pub fn cigar_string(&self) -> String {
        self.decomposition
            .iter()
            .flat_map(|d| d.cigar.clone())
            .map(|i| i.to_cigar_string())
            .collect::<Vec<String>>()
            .join("")
    }

    pub fn items(&self, seq: &[u8]) -> Vec<(Vec<u8>, Vec<CigarItem>, Vec<u8>)> {
        self.decomposition
            .iter()
            .map(|d| {
                (
                    self.motif_set.motifs[d.motif_idx].clone(),
                    d.cigar.clone(),
                    seq[d.start..d.end].to_vec(),
                )
            })
            .collect()
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

        for d in self.decomposition.iter() {
            let m = &self.motif_set.motifs[d.motif_idx];
            let mut qi = 0;
            let mut ri = d.start;
            for item in d.cigar.iter() {
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
    motif_alignment_score_cutoff: Option<i32>,
}

/// Representation of a motif alignment to a sequence.
/// Format: start (inclusive 0-based), end (inclusive 0-based), score, alignment table index
#[derive(Clone, Debug)]
pub struct MotifAlignmentInterval {
    start: usize,
    end: usize,
    score: i32,
    pub cigar: Vec<CigarItem>,
    motif_idx: usize,
}

fn backtrack_schedule(
    final_schedule: &mut Vec<MotifAlignmentInterval>,
    mut intervals: Vec<MotifAlignmentInterval>,
    m: &Vec<i32>,
    p: &Vec<usize>,
    j: usize,
) {
    if j == 0 {
        return;
    }
    if intervals[j - 1].score + m[p[j]] >= m[j - 1] {
        // avoid clone by moving intervals[j - 1]; we're done with everything from [j - 1] forward anyway
        // since we only ever decrement j (via p or via subtracting 1).
        let next_item = {
            let mut d = intervals.drain((j - 1)..);
            d.next().expect("backtract_schedule: next item must exist")
        };
        final_schedule.push(next_item);
        backtrack_schedule(final_schedule, intervals, m, p, p[j])
    } else {
        backtrack_schedule(final_schedule, intervals, m, p, j - 1)
    }
}

/// Implementation of known weighted interval scheduling algorithm to do the motif decomposition
/// See https://en.wikipedia.org/wiki/Interval_scheduling#Weighted
/// Returns the best (or one of the best) schedules + its score.
fn schedule(mut intervals: Vec<MotifAlignmentInterval>) -> (Vec<MotifAlignmentInterval>, i32) {
    // sort intervals globally by earliest to latest end index
    intervals.sort_by(|x, y| x.end.cmp(&y.end));

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
    fn to_cigar_string(&self) -> String {
        match self {
            Self::Ins(c) => format!("{}I", c),
            Self::Del(c) => format!("{}D", c),
            Self::Match(c) => format!("{}=", c),
            Self::Mismatch(c) => format!("{}X", c),
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

fn alignment_items_to_cigar(items: &[AlignmentItem]) -> Vec<CigarItem> {
    let mut cigar = Vec::new();
    let mut ii = items.iter();
    if let Some(first) = ii.next() {
        let mut current_op = first;
        let mut current_count: usize = 1;
        for op in ii {
            if op == current_op {
                current_count += 1;
            } else {
                cigar.push(current_op.to_cigar_item(current_count));
                current_op = op;
                current_count = 1;
            }
        }
        cigar.push(current_op.to_cigar_item(current_count));
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
            motif_set: Arc::new(motif_set),
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
    ) -> Option<(usize, usize, i32, Vec<CigarItem>)> {
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
            motif_alignment.push(if motif[row] == seq[col] {
                AlignmentItem::Match
            } else {
                AlignmentItem::Mismatch
            });
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

        //  2. determine intervals, cutting off low-scoring possibilities
        let intervals = self.compute_intervals(seq, &alignments, self.motif_alignment_score_cutoff, seq.len())?;

        //  3: use weighted interval scheduling algorithm https://en.wikipedia.org/wiki/Interval_scheduling#Weighted
        //     to find best sequence of motifs, with any 'idle' time being non-motif DNA in between motifs.
        let (decomposition, score) = schedule(intervals);

        Ok(MotifSequenceDecomposition {
            motif_set: self.motif_set.clone(),
            decomposition,
            score,
        })
    }

    pub fn decomp_to_str(&self, decomp: &MotifSequenceDecomposition) -> Result<Vec<&str>, Utf8Error> {
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
    #[case(b"CAG".to_vec(), vec!["CAG"], "CAG\n|||\nCAG")]
    #[case(b"CAAG".to_vec(), vec!["CAG"], "CA-G\n|| |\nCAAG")]
    #[case(
        b"CAGCAGCAGCAGCAGCAGCAGCAGCAG".to_vec(),
        vec!["CAG", "CAG", "CAG", "CAG", "CAG", "CAG", "CAG", "CAG", "CAG"],
        "CAGCAGCAGCAGCAGCAGCAGCAGCAG\n\
         |||||||||||||||||||||||||||\n\
         CAGCAGCAGCAGCAGCAGCAGCAGCAG",
    )]
    #[case(
        b"CCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCAG".to_vec(),
        vec!["CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG", "CCG",
             "CAG"],
        "CCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCAG\n\
         ||||||||||||||||||||||||||||||||||||||||||||||||\n\
         CCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCAG",
    )]
    #[case(
        b"CAGCAGCGGCAGCAAG".to_vec(),
        vec!["CAG", "CAG", "CGG", "CAG", "CAAG"],
        "CAGCAGCCGCAGCA-G\n\
         |||||||X|||||| |\n\
         CAGCAGCGGCAGCAAG",
    )]
    #[case(
        b"CAGCAGCAAGTTCAGCCGCCGCCCG".to_vec(),
        vec!["CAG", "CAG", "CAAG", "TT", "CAG", "CCG", "CCG", "CCCG"],
        "CAGCAGCA-G--CAGCCGCCGCC-G\n\
         |||||||| |  ||||||||||| |\n\
         CAGCAGCAAGTTCAGCCGCCGCCCG",
    )]
    fn test_decomposition(#[case] seq: Vec<u8>, #[case] expected_decomp: Vec<&str>, #[case] expected_align_str: &str) {
        let motif_set = MotifSet::new_from_strs(&vec!["CAG", "CCG"]);
        let decomposer = MotifSequenceDecomposer::new(motif_set, 5, -7, 4, Some(1));
        let res = decomposer.decompose(seq.as_slice()).unwrap();
        // assert_eq!(decomposer.decomp_to_str(&res).unwrap(), expected_decomp);
        assert_eq!(res.alignment_string(&seq), expected_align_str);
    }
}
