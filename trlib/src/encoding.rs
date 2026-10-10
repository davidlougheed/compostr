use rapidhash::RapidHashMap;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum EncodeError {
    #[error("alphabet exhausted")]
    AlphabetExhausted,
}

pub enum EncodeAlphabet {
    NewAlphabet(String),
}

fn encode_motif_sequence(motifs: &[&[u8]], alphabet: Vec<char>) -> Result<Vec<char>, EncodeError> {
    let mut res = Vec::with_capacity(motifs.len());
    let mut lookup: RapidHashMap<&[u8], char> = RapidHashMap::default();
    let mut alpha_ptr = 0usize;
    for &motif in motifs.iter() {
        if let Some(&val) = lookup.get(motif) {
            res.push(val);
        } else if alpha_ptr < alphabet.len() {
            let val = alphabet[alpha_ptr];
            lookup.insert(motif, val);
            alpha_ptr += 1;
            res.push(val);
        } else {
            return Err(EncodeError::AlphabetExhausted);
        }
    }
    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(vec!["CAG", "CAG", "CCG", "CAG"], "AB", "AABA")]
    #[case(vec!["CAG", "CCG", "CAG", "CCG"], "AB", "ABAB")]
    fn test_encoding(#[case] seq: Vec<&str>, #[case] alphabet: &str, #[case] end_res: &str) {
        let seq_u8: Vec<&[u8]> = seq.into_iter().map(|s| s.as_bytes()).collect();
        let res = encode_motif_sequence(&seq_u8, alphabet.chars().collect()).unwrap();
        assert_eq!(res, end_res.chars().collect::<Vec<char>>());
    }
}
