pub struct MotifSet {
    pub motifs: Vec<Vec<u8>>,
}

impl MotifSet {
    pub fn new(motifs: Vec<Vec<u8>>) -> Self {
        MotifSet { motifs }
    }

    pub fn new_from_strs(motifs: &[&str]) -> Self {
        MotifSet::new(motifs.iter().map(|&m| m.bytes().collect()).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_motif_set() {
        MotifSet::new(vec![b"CAG".to_vec(), b"CCG".to_vec()]);
        MotifSet::new_from_strs(&vec!["CAG", "CCG"]);
    }
}
