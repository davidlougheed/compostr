use std::fs::File;
use std::io::{BufReader, Error as IoError, Write};

use clap::Parser;
use noodles_core::Position;
use noodles_fasta as fasta;
use serde::Serialize;
use serde_json;
use thiserror::Error;
use trlib::decomposition::{MotifSequenceDecomposer, MotifSequenceDecomposition};
use trlib::motif::MotifSet;

use crate::CompostrCliError::FastaIoError;

#[derive(Debug, Error)]
pub enum CompostrCliError {
    #[error("FASTA IO error: {0}")]
    FastaIoError(IoError),
}

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[arg(short, long, num_args=0..)]
    motifs: Vec<String>,

    #[arg(long)]
    out: String,

    fasta: String,
}

#[derive(Serialize)]
struct OutputResult {
    name: String,
    sequence: String,
    result: MotifSequenceDecomposition,
}

fn process_fasta(cli: &Cli) -> Result<Vec<OutputResult>, CompostrCliError> {
    let motifs: Vec<&str> = cli.motifs.iter().map(|i| i.as_str()).collect();
    let motif_set = MotifSet::new_from_strs(&motifs);
    // TODO: parameterize
    let decomposer = MotifSequenceDecomposer::new(motif_set, 4, -4, 5, Some(-1)).unwrap();

    let mut reader = File::open(&cli.fasta)
        .map(BufReader::new)
        .map(fasta::io::Reader::new)
        .expect("could not open FASTA");

    let start = Position::try_from(1).unwrap();

    let mut outputs = Vec::new();

    for result in reader.records() {
        let record = result.map_err(|e| FastaIoError(e))?;
        let seq_opt = record.sequence().get(start..);
        if let Some(seq) = seq_opt {
            let res = decomposer.decompose(seq).unwrap();
            outputs.push(OutputResult {
                name: str::from_utf8(record.name()).unwrap().to_string(),
                sequence: str::from_utf8(seq).unwrap().to_string(),
                result: res,
            });
        }
    }

    Ok(outputs)
}


fn main() -> Result<(), CompostrCliError> {
    let cli = Cli::parse();
    let outputs = process_fasta(&cli)?;

    let json_res = serde_json::to_string_pretty(&outputs).unwrap();
    let mut writer = File::create(&cli.out).expect("could not open outfile");
    writer.write(json_res.as_bytes()).expect("could not write to outfile");

    Ok(())
}
