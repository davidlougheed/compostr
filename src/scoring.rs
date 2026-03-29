use std::borrow::Borrow;
use std::collections::HashMap;

use once_cell::sync::Lazy;
use parasail_rs::prelude::{Error, Matrix};

fn make_iupac_lookup(t_base: u8) -> HashMap<u8, Vec<u8>> {
    let mut lookup = HashMap::new();
    lookup.insert(b'A', vec![b'A']);
    lookup.insert(b'C', vec![b'C']);
    lookup.insert(b'G', vec![b'G']);
    lookup.insert(t_base, vec![t_base]);
    lookup.insert(b'M', vec![b'A', b'C']);
    lookup.insert(b'R', vec![b'A', b'G']);
    lookup.insert(b'W', vec![b'A', t_base]);
    lookup.insert(b'S', vec![b'C', b'G']);
    lookup.insert(b'Y', vec![b'C', t_base]);
    lookup.insert(b'K', vec![b'G', t_base]);
    lookup.insert(b'V', vec![b'A', b'C', b'G']);
    lookup.insert(b'H', vec![b'A', b'C', t_base]);
    lookup.insert(b'D', vec![b'A', b'G', t_base]);
    lookup.insert(b'B', vec![b'C', b'G', t_base]);
    lookup.insert(b'N', vec![b'A', b'C', b'G', t_base]);
    lookup
}

const IUPAC_CODES_DNA: &[u8; 15] = b"ACGTRYSWKMBDHVN";
static IUPAC_CODES_DNA_REVERSE: Lazy<HashMap<u8, usize>> =
    Lazy::new(|| HashMap::from_iter(IUPAC_CODES_DNA.iter().enumerate().map(|(i, &b)| (b, i))));

static IUPAC_CODE_DNA_LOOKUP: Lazy<HashMap<u8, Vec<u8>>> = Lazy::new(|| make_iupac_lookup(b'T'));

#[derive(Debug)]
pub enum ScoringMatrixError {
    ParasailError(Error),
}

/// Right now essentially just a wrapper for parasail-rs's Matrix struct so that we can define the API.
pub struct ScoringMatrix {
    pub matrix: Matrix,
}

impl ScoringMatrix {
    pub fn new_iupac_dna(match_score: i32, mismatch_score: i32) -> Result<Self, ScoringMatrixError> {
        let mut matrix = Matrix::create(IUPAC_CODES_DNA, match_score, mismatch_score)
            .map_err(|e| ScoringMatrixError::ParasailError(e))?;

        for (code, code_matches) in IUPAC_CODE_DNA_LOOKUP.borrow().iter() {
            for cm in code_matches.iter() {
                matrix
                    .set_value(
                        IUPAC_CODES_DNA_REVERSE[code] as i32,
                        IUPAC_CODES_DNA_REVERSE[cm] as i32,
                        match_score,
                    )
                    .map_err(|e| ScoringMatrixError::ParasailError(e))?;
                matrix
                    .set_value(
                        IUPAC_CODES_DNA_REVERSE[cm] as i32,
                        IUPAC_CODES_DNA_REVERSE[code] as i32,
                        match_score,
                    )
                    .map_err(|e| ScoringMatrixError::ParasailError(e))?;
            }
        }

        Ok(ScoringMatrix { matrix })
    }
}
