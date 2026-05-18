use trlib::decomposition::MotifSequenceDecomposer;
use trlib::motif::MotifSet;

fn main() {
    let seq = b"CAGCAGCAGCAGCAG".to_vec();
    let motif_set = MotifSet::new_from_strs(&vec!["CAG"]);
    let decomposer = MotifSequenceDecomposer::new(motif_set, 4, -4, 5, Some(-1)).unwrap();
    let res = decomposer.decompose(seq.as_slice()).unwrap();
    eprintln!("{}", res.alignment_string(&seq));

    let seq2 = b"TTTTTATTTTTATTTTTTATTTTTCTT".to_vec();
    let motif_set2 = MotifSet::new_from_strs(&vec!["TTTTTAT"]);
    let decomposer2 = MotifSequenceDecomposer::new(motif_set2, 4, -4, 5, Some(-1)).unwrap();
    let res2 = decomposer2.decompose(seq2.as_slice()).unwrap();
    eprintln!("{}", res2.alignment_string(&seq2));
}
