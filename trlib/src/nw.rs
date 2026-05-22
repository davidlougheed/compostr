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
    m_gap_len: Array2<i32>,
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

    pub fn align_motif_to_seq(self, motif: &[u8], seq: &[u8]) -> Alignment {
        let seq_len = seq.len();
        let motif_len = motif.len();
        let mut m_score: Array2<i32> = Array::zeros((seq_len + 1, motif_len + 1));

        let mut m_trace: Array2<TraceItem> = Array::from_elem((seq.len() + 1, motif.len() + 1), TraceItem::Unset);
        m_trace[[0, 0]] = TraceItem::Done;

        let mut m_gap_len: Array2<i32> = Array::zeros((seq_len + 1, motif_len + 1));

        for i in 1..seq_len {
            let mut x = i;
            let mut y = 1;
            //iterate the diagonal
            while y <= motif_len && x > 0 {
                //up
                let up = match m_trace[[x,y-1]] {
                    TraceItem::Up | TraceItem::Left => m_score[[x,y-1]] + self.gap_extend,
                    _ => m_score[[x,y-1]] + self.gap_open
                };
                //left
                let left = match m_trace[[x-1,y]] {
                    TraceItem::Up | TraceItem::Left => m_score[[x-1,y]] + self.gap_extend,
                    _ => m_score[[x-1,y]] + self.gap_open
                };
                //(mis)match - upleft
                let upleft = if motif[y] == seq[x] {
                    m_score[[x-1,y-1]] + self.match_score
                } else {
                    m_score[[x-1,y-1]] + self.mismatch_score
                };

                if up > upleft && up > left {
                    m_score[[x,y]] = up;
                    m_trace[[x,y]] = TraceItem::Up;
                    m_gap_len[[x,y]] = m_gap_len[[x,y-1]] + 1;
                } else if left > upleft && left > up {
                    m_score[[x,y]] = left;
                    m_trace[[x,y]] = TraceItem::Left;
                    m_gap_len[[x,y]] = m_gap_len[[x-1,y]] + 1;
                } else {
                    m_score[[x,y]] = upleft;
                    m_trace[[x,y]] = TraceItem::Diag;
                }
                x -= 1;
                y += 1;
            }
        }

        Alignment { m_score, m_trace, m_gap_len}
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
}
