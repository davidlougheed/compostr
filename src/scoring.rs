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

fn rev_from_alphabet(alphabet: &[u8]) -> HashMap<u8, usize> {
    HashMap::from_iter(alphabet.iter().enumerate().map(|(i, &b)| (b, i)))
}

const IUPAC_DNA_ALPHABET: &[u8; 15] = b"ACGTRYSWKMBDHVN";
static IUPAC_DNA_ALPHABET_REVERSE: Lazy<HashMap<u8, usize>> = Lazy::new(|| rev_from_alphabet(IUPAC_DNA_ALPHABET));

const IUPAC_RNA_ALPHABET: &[u8; 15] = b"ACGURYSWKMBDHVN";
static IUPAC_RNA_ALPHABET_REVERSE: Lazy<HashMap<u8, usize>> = Lazy::new(|| rev_from_alphabet(IUPAC_RNA_ALPHABET));

/// Lookup from IUPAC code to canonical DNA base
static IUPAC_CODE_DNA_LOOKUP: Lazy<HashMap<u8, Vec<u8>>> = Lazy::new(|| make_iupac_lookup(b'T'));
/// Lookup from IUPAC code to canonical RNA base
static IUPAC_CODE_RNA_LOOKUP: Lazy<HashMap<u8, Vec<u8>>> = Lazy::new(|| make_iupac_lookup(b'U'));

#[derive(Debug)]
pub enum ScoringMatrixError {
    ParasailError(Error),
}

/// Right now essentially just a wrapper for parasail-rs's Matrix struct so that we can define the API.
pub struct ScoringMatrix {
    pub matrix: Matrix,
}

fn make_parasail_matrix(is_rna: bool, match_score: i32, mismatch_score: i32) -> Result<Matrix, Error> {
    let alphabet = if is_rna { IUPAC_RNA_ALPHABET } else { IUPAC_DNA_ALPHABET };
    let mut matrix = Matrix::create(alphabet, match_score, mismatch_score)?;

    let code_lookup = if is_rna {
        IUPAC_CODE_RNA_LOOKUP.borrow()
    } else {
        IUPAC_CODE_DNA_LOOKUP.borrow()
    };

    let rev = if is_rna {
        IUPAC_RNA_ALPHABET_REVERSE.borrow()
    } else {
        IUPAC_DNA_ALPHABET_REVERSE.borrow()
    };

    for (code, code_matches) in code_lookup.iter() {
        for cm in code_matches.iter() {
            matrix.set_value(rev[code] as i32, rev[cm] as i32, match_score)?;
            matrix.set_value(rev[cm] as i32, rev[code] as i32, match_score)?;
        }
    }

    Ok(matrix)
}

impl ScoringMatrix {
    pub fn new_iupac_dna(match_score: i32, mismatch_score: i32) -> Result<Self, ScoringMatrixError> {
        Ok(ScoringMatrix {
            matrix: make_parasail_matrix(false, match_score, mismatch_score)
                .map_err(|e| ScoringMatrixError::ParasailError(e))?,
        })
    }

    pub fn new_iupac_rna(match_score: i32, mismatch_score: i32) -> Result<Self, ScoringMatrixError> {
        Ok(ScoringMatrix {
            matrix: make_parasail_matrix(true, match_score, mismatch_score)
                .map_err(|e| ScoringMatrixError::ParasailError(e))?,
        })
    }
}
