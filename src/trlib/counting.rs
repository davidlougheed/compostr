use parasail_rs::prelude::{Aligner, Alignment, Error};

use crate::motif::MotifSet;

/// TODO
pub struct MotifSequenceDecomposition<'m> {
    decomposition: Vec<&'m [u8]>,
}

/// TODO
pub struct MotifSequenceDecomposer {
    motif_set: MotifSet,
}

/// TODO
fn schedule<'m>(alignments: &Vec<(&'m [u8], Alignment)>) -> Vec<&'m [u8]> {
    // re-use known weighted interval scheduling algorithm to do the motif decomposition
    // TODO
    Vec::new()
}

impl MotifSequenceDecomposer {
    /// TODO
    pub fn decompose(&self, seq: &[u8]) -> Result<MotifSequenceDecomposition, Error> {
        // rough algorithm outline, 2 parts:

        //  1: align all motifs (ends-free) to sequence to get alignment score matrix
        let aligner = Aligner::new().semi_global().striped().build();
        let mut alignments: Vec<(&[u8], Alignment)> = Vec::with_capacity(self.motif_set.motifs.len());
        for m in self.motif_set.motifs.iter() {
            alignments.push((m, aligner.align(Some(m), seq)?));
        }

        //  2. determine intervals using some kind of heuristic so we don't have an absurd number?
        //     or just use last row(?) of the matrix as the score + figure out the interval... + do a little trimming

        //  3: use weighted interval scheduling algorithm https://en.wikipedia.org/wiki/Interval_scheduling#Weighted
        //     to find best sequence of motifs, with any 'idle' time being non-motif DNA in between motifs.
        let decomposition = schedule(&alignments);

        Ok(MotifSequenceDecomposition { decomposition })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
}
