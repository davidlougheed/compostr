use std::sync::{Arc, Mutex};

use rapidhash::{HashMapExt, RapidHashMap};
use thiserror::Error;

use crate::decomposition::MotifSequenceDecomposition;

#[derive(Error, PartialEq, Debug)]
pub enum EncodeError {
    #[error("alphabet exhausted")]
    AlphabetExhausted,
    #[error("could not secure mutex lock")]
    MutexLockError,
}

#[derive(Debug)]
pub struct EncodeAlphabet {
    pub alphabet: Vec<char>,
    ptr: usize,
    // could maybe make keys slices if we could figure out lifetimes
    lookup: RapidHashMap<Vec<u8>, char>,
}

impl EncodeAlphabet {
    /// Returns an Option of the char representing the specified motif.
    pub fn get_motif_code(&self, motif: &[u8]) -> Option<char> {
        self.lookup.get(motif).cloned()
    }

    /// Encodes the specified motif, advancing the pointer and allocating a new code from the alphabet if needed.
    /// If the alphabet is exhausted, an error will be returned.
    pub fn encode_motif(&mut self, motif: &[u8]) -> Result<char, EncodeError> {
        if let Some(&val) = self.lookup.get(motif) {
            Ok(val)
        } else if self.ptr < self.alphabet.len() {
            let val = self.alphabet[self.ptr];
            self.lookup.insert(motif.to_vec(), val);
            self.ptr += 1;
            Ok(val)
        } else {
            Err(EncodeError::AlphabetExhausted)
        }
    }
}

pub enum EncodeAlphabetParam {
    NewAlphabet(Vec<char>),
    ExistingAlphabet(Arc<Mutex<EncodeAlphabet>>),
}

pub struct Encoder {
    alphabet: EncodeAlphabetParam,
}

impl Encoder {
    pub fn new(alphabet: EncodeAlphabetParam) -> Self {
        Encoder { alphabet }
    }

    pub fn new_from_str(alphabet_str: &str) -> Self {
        let alphabet = EncodeAlphabetParam::NewAlphabet(alphabet_str.chars().collect());
        Encoder { alphabet }
    }

    pub fn get_alphabet_chars(&self) -> Result<Vec<char>, EncodeError> {
        match &self.alphabet {
            EncodeAlphabetParam::NewAlphabet(a) => Ok(a.clone()),
            EncodeAlphabetParam::ExistingAlphabet(a) => {
                let al = a.lock().map_err(|_| EncodeError::MutexLockError)?;
                Ok(al.alphabet.clone())
            }
        }
    }

    pub fn encode(&mut self, motifs: &[&[u8]]) -> Result<(Vec<char>, Arc<Mutex<EncodeAlphabet>>), EncodeError> {
        let a = match &self.alphabet {
            EncodeAlphabetParam::NewAlphabet(a) => Arc::new(Mutex::new(EncodeAlphabet {
                alphabet: a.clone(),
                ptr: 0,
                lookup: RapidHashMap::new(),
            })),
            EncodeAlphabetParam::ExistingAlphabet(a) => a.clone(),
        };
        let mut al = a.lock().map_err(|_| EncodeError::MutexLockError)?;

        let res = motifs.iter().map(|&motif| al.encode_motif(motif)).collect::<Result<Vec<_>, _>>()?;

        // Needs to occur after all encode results are collected, since that mutates the alphabet.
        self.alphabet = EncodeAlphabetParam::ExistingAlphabet(a.clone());

        Ok((res, a.clone()))
    }

    pub fn encode_vec_u8(&mut self, motifs: Vec<Vec<u8>>) -> Result<(Vec<char>, Arc<Mutex<EncodeAlphabet>>), EncodeError> {
        let motif_refs = motifs.iter().map(|m| m.as_slice()).collect::<Vec<&[u8]>>();
        self.encode(&motif_refs)
    }

    pub fn encode_decomposition(
        &mut self,
        seq: &[u8],
        d: &MotifSequenceDecomposition,
        with_interruptions: bool,
    ) -> Result<(Vec<char>, Arc<Mutex<EncodeAlphabet>>), EncodeError> {
        let si = d.sequence_items(seq, with_interruptions);
        self.encode(&si)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(vec!["CAG", "CAG", "CCG", "CAG"], "AB", "AABA")]
    #[case(vec!["CAG", "CCG", "CAG", "CCG"], "AB", "ABAB")]
    fn test_encoding(#[case] seq: Vec<&str>, #[case] alphabet: &str, #[case] end_res: &str) {
        let alphabet_chars: Vec<char> = alphabet.chars().collect();
        let mut encoder = Encoder {
            alphabet: EncodeAlphabetParam::NewAlphabet(alphabet_chars),
        };
        let seq_u8: Vec<&[u8]> = seq.into_iter().map(|s| s.as_bytes()).collect();
        let res = encoder.encode(&seq_u8).expect("should work");
        assert_eq!(res.0, end_res.chars().collect::<Vec<char>>());

        let res2 = encoder.encode(&seq_u8).expect("should work again");
        assert_eq!(res.0, res2.0);
    }

    #[rstest]
    #[case(vec!["CAG", "CAG", "CCG", "CGG"], "AB")]
    fn test_encoding_exhausted(#[case] seq: Vec<&str>, #[case] alphabet: &str) {
        let mut encoder = Encoder {
            alphabet: EncodeAlphabetParam::NewAlphabet(alphabet.chars().collect()),
        };
        let seq_u8: Vec<&[u8]> = seq.into_iter().map(|s| s.as_bytes()).collect();
        let err = encoder.encode(&seq_u8).expect_err("should be exhausted");
        assert_eq!(err, EncodeError::AlphabetExhausted)
    }

    #[rstest]
    #[case(vec!["CAG", "CAG", "CCG", "CCG"], "AB", "AABB", vec!["CAG", "CCG", "CCG", "CCG"], "ABBB")]
    fn test_encoding_reuse(
        #[case] seq: Vec<&str>,
        #[case] alphabet: &str,
        #[case] end_res_1: &str,
        #[case] seq_2: Vec<&str>,
        #[case] end_res_2: &str,
    ) {
        let mut encoder = Encoder {
            alphabet: EncodeAlphabetParam::NewAlphabet(alphabet.chars().collect()),
        };
        let seq_u8: Vec<&[u8]> = seq.into_iter().map(|s| s.as_bytes()).collect();

        let res1 = encoder.encode(&seq_u8).expect("res1 should work");
        assert_eq!(res1.0, end_res_1.chars().collect::<Vec<char>>());

        let mut encoder_2 = Encoder {
            alphabet: EncodeAlphabetParam::ExistingAlphabet(res1.1),
        };

        let seq_2_u8: Vec<&[u8]> = seq_2.into_iter().map(|s| s.as_bytes()).collect();
        let res2 = encoder_2.encode(&seq_2_u8).expect("res2 should work");
        assert_eq!(res2.0, end_res_2.chars().collect::<Vec<char>>());
    }
}
