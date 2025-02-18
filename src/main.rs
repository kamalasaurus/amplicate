use std::env::args;
use std::fs::File;
use paraseq::{
    fasta,
    fastx::Record,
    parallel::{ParallelProcessor, ParallelReader, ProcessError}
};

#[derive(Clone, Default)]
struct MyProcessor {
    // Add fields here
}

impl ParallelProcessor for MyProcessor {
    fn process_record<R: Record>(&mut self, record: R) -> Result<(), ProcessError> {
        // Process record in parallel
        println!("Processing record: {}", record.id_str());
        println!("Sequence: {}", record.seq_str());
        println!("Length: {}", record.seq().len());
        Ok(())
    }
}

fn main() -> Result<(), ProcessError> {
    let input = args().nth(1).expect("No input file provided.");
    let file = File::open(input)?;
    let reader = fasta::Reader::new(file);
    let processor = MyProcessor::default();
    let num_threads = 8;
    reader.process_parallel(processor, num_threads)?;
    Ok(())
}