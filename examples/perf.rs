use std::time::SystemTime;

use trlib::trlib::decomposition::MotifSequenceDecomposer;
use trlib::trlib::motif::MotifSet;

fn main() {
    let motif_set = MotifSet::new_from_strs(&vec!["CAG", "CCG"]);
    let seq = b"CCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCAG".to_vec();
    let d = MotifSequenceDecomposer::new(motif_set, 5, -7, 4, Some(1));

    let t = SystemTime::now();
    for _ in 1..100000 {
        d.decompose(&seq).unwrap();
    }
    println!("{}µs", SystemTime::now().duration_since(t).unwrap().as_micros());
}
