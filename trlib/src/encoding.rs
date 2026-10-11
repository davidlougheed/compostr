use std::sync::{Arc, Mutex};

use rapidhash::{HashMapExt, RapidHashMap};
use thiserror::Error;

use crate::decomposition::MotifSequenceDecomposition;

#[derive(Error, Debug)]
pub enum EncodeError {
    #[error("alphabet exhausted")]
    AlphabetExhausted,
    #[error("could not secure mutex lock")]
    MutexLockError,
}

pub struct EncodeAlphabet<'a> {
    alphabet: &'a [char],
    ptr: usize,
    lookup: RapidHashMap<&'a [u8], char>,
}

impl<'a> EncodeAlphabet<'a> {
    pub fn encode_motif(&mut self, motif: &'a [u8]) -> Result<char, EncodeError> {
        if let Some(&val) = self.lookup.get(motif) {
            Ok(val)
        } else if self.ptr < self.alphabet.len() {
            let val = self.alphabet[self.ptr];
            self.lookup.insert(motif, val);
            self.ptr += 1;
            Ok(val)
        } else {
            Err(EncodeError::AlphabetExhausted)
        }
    }
}

pub enum EncodeAlphabetParam<'a> {
    NewAlphabet(&'a [char]),
    ExistingAlphabet(Arc<Mutex<EncodeAlphabet<'a>>>),
}

pub struct Encoder<'a> {
    alphabet: EncodeAlphabetParam<'a>,
}

impl<'a> Encoder<'a> {
    pub fn encode(
        &'a mut self,
        motifs: &[&'a [u8]],
    ) -> Result<(Vec<char>, Arc<Mutex<EncodeAlphabet<'a>>>), EncodeError> {
        let mut res = Vec::with_capacity(motifs.len());
        let a = match &self.alphabet {
            EncodeAlphabetParam::NewAlphabet(a) => Arc::new(Mutex::new(EncodeAlphabet {
                alphabet: a,
                ptr: 0,
                lookup: RapidHashMap::new(),
            })),
            EncodeAlphabetParam::ExistingAlphabet(a) => a.clone(),
        };
        let mut al = a.lock().map_err(|_| EncodeError::MutexLockError)?;

        for &motif in motifs.iter() {
            res.push(al.encode_motif(motif)?);
        }

        self.alphabet = EncodeAlphabetParam::ExistingAlphabet(a.clone());

        Ok((res, a.clone()))
    }

    pub fn encode_decomposition(
        &'a mut self,
        seq: &'a [u8],
        d: &'a MotifSequenceDecomposition,
        with_interruptions: bool,
    ) -> Result<(Vec<char>, Arc<Mutex<EncodeAlphabet<'a>>>), EncodeError> {
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
            alphabet: EncodeAlphabetParam::NewAlphabet(&alphabet_chars),
        };
        let seq_u8: Vec<&[u8]> = seq.into_iter().map(|s| s.as_bytes()).collect();
        let res = encoder.encode(&seq_u8).expect("should work");
        assert_eq!(res.0, end_res.chars().collect::<Vec<char>>());
    }
}
