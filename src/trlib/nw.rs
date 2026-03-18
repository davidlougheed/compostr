use ndarray::{Array, Array2};

#[derive(Clone)]
pub enum TraceItem {
    Unset,
    Done,
    Up,
    Left,
    Diag,
}

pub struct Alignment {
    m_score: Array2<i32>,
    m_trace: Array2<TraceItem>,
}

pub struct Aligner {
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
}

impl Aligner {
    pub fn new(match_score: i32, mismatch_score: i32, gap_open: i32, gap_extend: i32) -> Self {
        Aligner {
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
        }
    }

    pub fn align_motif_to_seq(motif: &[u8], seq: &[u8]) -> Alignment {
        let mut m_score: Array2<i32> = Array::zeros((seq.len() + 1, motif.len() + 1));
        // TODO: first col with gaps

        let mut m_trace: Array2<TraceItem> = Array::from_elem((seq.len() + 1, motif.len() + 1), TraceItem::Unset);
        m_trace[[0, 0]] = TraceItem::Done;

        // for i in 0..seq.len() {
        //     for
        // }

        Alignment { m_score, m_trace }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
}
