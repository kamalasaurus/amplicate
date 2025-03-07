use std::env::args;
use std::fs::File;
use std::collections::HashMap;
use paraseq::{
    fasta,
    fastx::Record,
    parallel::{ParallelProcessor, ParallelReader, ProcessError}
};
use num_cpus;

// Helper function: finds needle in haystack
// using case insensitive matching, with 'N' as wildcard.
fn find_with_wildcard(haystack: &str, needle: &str) -> Option<usize> {
    let haystack = haystack.as_bytes();
    let needle = needle.as_bytes();
    if needle.len() > haystack.len() { return None; }
    for i in 0..=haystack.len() - needle.len() {
        let mut found = true;
        for j in 0..needle.len() {
            let h = haystack[i + j].to_ascii_uppercase();
            let n = needle[j].to_ascii_uppercase();
            if h != n && h != b'N' && n != b'N' {
                found = false;
                break;
            }
        }
        if found { return Some(i); }
    }
    None
}

// Update MyProcessor to store shared primers and an optional prefix length.
#[derive(Clone)]
struct MyProcessor {
    primers: HashMap<String, String>, // key: primer id, value: sequence
    prefix_len: Option<usize>,
}

impl Default for MyProcessor {
    fn default() -> Self {
        Self { 
            primers: HashMap::new(),
            prefix_len: None,
        }
    }
}

impl ParallelProcessor for MyProcessor {
    fn process_record<R: Record>(&mut self, record: R) -> Result<(), ProcessError> {
        let record_id = record.id_str();
        let seq = record.seq_str();
        let primer_vec: Vec<(&String, &String)> = self.primers.iter().collect();
        for i in 0..primer_vec.len() {
            for j in (i + 1)..primer_vec.len() {
                let (id1, p1) = primer_vec[i];
                let (id2, p2) = primer_vec[j];
                if let Some(n) = self.prefix_len {
                    if id1.len() < n || id2.len() < n || &id1[..n] != &id2[..n] {
                        continue;
                    }
                }
                if let Some(pos1) = find_with_wildcard(seq, p1) {
                    if let Some(pos2) = find_with_wildcard(seq, p2) {
                        if pos1 <= pos2 {
                            // Include p1 and p2 in the amplicon.
                            let amplicon = &seq[pos1..(pos2 + p2.len())];
                            println!("{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                                record_id, id1, p1, id2, p2, amplicon, amplicon.len());
                        } else {
                            let amplicon = &seq[pos2..(pos1 + p1.len())];
                            println!("{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                                record_id, id2, p2, id1, p1, amplicon, amplicon.len());
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

fn main() -> Result<(), ProcessError> {
    let input_genome = args().nth(1).expect("No input genome provided.");
    let input_primers = args().nth(2).expect("No input primers provided.");
    // Optional flag: if provided, parse the prefix length N.
    let prefix_len = args().nth(3).and_then(|s| s.parse::<usize>().ok());
    
    // Load primers (shared for all threads) into a HashMap.
    let primers_file = File::open(input_primers)?;
    let mut primers_reader = fasta::Reader::new(primers_file);
    let mut primers: HashMap<String, String> = HashMap::new();
    let mut primer_set = fasta::RecordSet::new(5096); // buffer up to 5096 primers

    while primer_set.fill(&mut primers_reader)? {
        for record in primer_set.iter() {
            let record = record?;
            primers.insert(record.id_str().to_string(), record.seq_str().to_string());
        }
    }

    let genome_file = File::open(input_genome)?;
    let genome_reader = fasta::Reader::new(genome_file);
    let mut processor = MyProcessor { primers, prefix_len: None };
    // Update processor with prefix_len if available.
    if let Some(n) = prefix_len {
        processor.prefix_len = Some(n);
    }
    
    // Set num_threads to max(num_cpus - 2, 1)
    let num_threads = (num_cpus::get()).saturating_sub(2).max(1);
    
    genome_reader.process_parallel(processor, num_threads)?;
    Ok(())
}