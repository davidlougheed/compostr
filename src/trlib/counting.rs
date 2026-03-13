use std::str::{self, Utf8Error};
use parasail_rs::prelude::{Aligner, Alignment, Error};

use crate::motif::MotifSet;

/// TODO
pub struct MotifSequenceDecomposition {
    pub decomposition: Vec<Vec<u8>>,
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

/// TODO
pub struct MotifSequenceDecomposer {
    pub motif_set: MotifSet,
}

/// TODO
fn schedule(alignments: &Vec<(&[u8], Alignment)>) -> Vec<Vec<u8>> {
    // re-use known weighted interval scheduling algorithm to do the motif decomposition
    // https://en.wikipedia.org/wiki/Interval_scheduling#Weighted

    // TODO

    // build up a decomposition of motifs, our "schedule"
    let decomposition = Vec::new();

    // TODO

    decomposition
}

impl MotifSequenceDecomposer {
    pub fn new(motif_set: MotifSet) -> Self {
        MotifSequenceDecomposer { motif_set }
    }

    /// TODO
    pub fn decompose<'m>(&'m self, seq: &[u8]) -> Result<MotifSequenceDecomposition, Error> {
        // rough algorithm outline, 2 parts:

        //  1: align all motifs (ends-free) to sequence to get alignment score matrix
        let aligner = Aligner::new().semi_global().striped().use_table().build();
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

    #[test]
    fn test_decomposition() {
        let motif_set = MotifSet::new_from_strs(&vec!["CAG", "CCG"]);
        let decomposer = MotifSequenceDecomposer::new(motif_set);

        let res1 = decomposer.decompose(b"CAGCAGCAAGTTCAGCCGCCGCCCG").unwrap();
        assert_eq!(
            res1.decomposition_strs().unwrap(),
            vec!["CAG", "CAG", "CAAG", "T", "T", "CAG", "CCG", "CCG", "CCCG"]
        );
    }
}
