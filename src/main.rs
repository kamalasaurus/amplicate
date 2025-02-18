use std::env::args;
use bio::io::fasta;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = args().nth(1).unwrap();
    let reader = fasta::Reader::from_file(input).unwrap();
    reader
        .records()
        .try_for_each(|record| -> Result<(), Box<dyn std::error::Error>> {
            let record = record.map_err(|e| Box::new(e) as Box<dyn std::error::Error>)?;
            println!("Record ID: {}", record.id());
            println!("Description: {:?}", record.desc());
            println!("Sequence: {}", String::from_utf8_lossy(record.seq()));
            println!("Length: {}", record.seq().len());
            Ok(())
        })?;

    Ok(())
}