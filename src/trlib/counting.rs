use crate::motif::MotifSet;

///
pub struct MotifSequenceDecomposition {

}

pub struct MotifSequenceDecomposer {
    motif_set: MotifSet,
}

impl MotifSequenceDecomposer {
    fn decompose(&self, seq: &str) {
        // rough algorithm outline, 2 parts:
        //  1: align all motifs (ends-free) to sequence to get alignment score matrix
        //  2: use weighted interval scheduling algorithm https://en.wikipedia.org/wiki/Interval_scheduling#Weighted
        //     to find best sequence of motifs, with any 'idle' time being non-motif DNA in between motifs.
    }
}

#[cfg(test)]
mod tests {
    use super::*;
}
