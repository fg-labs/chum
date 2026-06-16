//! BLASTn subprocess wrapper, format-6 parsing, and hit summarization.
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::thread;

use anyhow::{Context, Result, bail};
use serde::Deserialize;

use crate::intervals::Bait;

/// BLASTn format-6 output record (12 standard columns).
#[derive(Debug, Clone, Deserialize)]
pub struct BlastHitFormat6 {
    /// Query sequence identifier.
    pub qseqid: String,
    /// Subject (database) sequence identifier.
    pub sseqid: String,
    /// Percentage of identical positions (0–100).
    pub pident: f64,
    /// Alignment length in bases.
    pub length: u32,
    /// Number of mismatched positions.
    pub mismatch: u32,
    /// Number of gap-open events.
    pub gapopen: u32,
    /// Start position of the alignment on the query (1-based).
    pub qstart: u32,
    /// End position of the alignment on the query (1-based).
    pub qend: u32,
    /// Start position of the alignment on the subject (1-based).
    pub sstart: u32,
    /// End position of the alignment on the subject (1-based).
    pub send: u32,
    /// Expect value (e-value) of the alignment.
    pub evalue: f64,
    /// Bit score of the alignment.
    pub bitscore: f64,
}

/// Configuration for a BLASTn run.
#[derive(Debug, Clone)]
pub struct BlastRunner {
    /// BLAST database name.
    pub db: String,
    /// Directory containing the database files.
    pub db_path: Option<PathBuf>,
    /// Whether to enable DUST low-complexity masking.
    pub dust: bool,
    /// Number of CPU threads passed to blastn via `-num_threads`.
    pub num_threads: usize,
}

impl BlastRunner {
    /// Create a new [`BlastRunner`].
    pub fn new(
        db: impl Into<String>,
        db_path: Option<PathBuf>,
        dust: bool,
        num_threads: usize,
    ) -> Self {
        BlastRunner {
            db: db.into(),
            db_path,
            dust,
            num_threads,
        }
    }

    /// Align a batch of bait sequences against the database.
    ///
    /// Each sequence is written to the query FASTA under its zero-based batch index
    /// (`>{i}`) rather than its bait name. Using the index as the identifier avoids
    /// collisions between baits that share a name (or have no name at all) and keeps
    /// every query uniquely addressable regardless of name or sequence content; hits are
    /// routed back to the originating bait by parsing `qseqid` as that index. The query
    /// FASTA is streamed to `blastn` via stdin and hits are read from stdout, avoiding
    /// temporary files entirely.
    ///
    /// Returns one `Vec<BlastHitFormat6>` per input bait, in the same order.
    pub fn align_batch(&self, baits: &[Bait]) -> Result<Vec<Vec<BlastHitFormat6>>> {
        if baits.is_empty() {
            return Ok(vec![]);
        }

        let mut query = Vec::new();
        for (i, bait) in baits.iter().enumerate() {
            if let Some(seq) = &bait.sequence {
                writeln!(query, ">{i}\n{seq}")?;
            } else {
                bail!("Bait '{}' has no sequence! Set one.", bait.name);
            }
        }

        // Spawn blastn: read query from stdin, write hits to stdout.
        let mut cmd = Command::new("blastn");
        cmd.arg("-query")
            .arg("-")
            .arg("-out")
            .arg("-")
            .arg("-outfmt")
            .arg("6")
            .arg("-db")
            .arg(&self.db)
            .arg("-dust")
            .arg(if self.dust { "yes" } else { "no" })
            .arg("-num_threads")
            .arg(self.num_threads.to_string())
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped());

        if let Some(db_path) = &self.db_path {
            cmd.env("BLASTDB", db_path);
        }

        let mut child = cmd
            .spawn()
            .with_context(|| "Failed to execute `blastn`; ensure it is on $PATH")?;

        // Write FASTA to stdin in a background thread to avoid deadlock when
        // blastn's stdout pipe buffer fills before we finish writing the query.
        let mut stdin = child.stdin.take().unwrap();
        let writer = thread::spawn(move || stdin.write_all(&query));

        let output = child
            .wait_with_output()
            .with_context(|| "Failed to wait for `blastn`")?;
        // Capture the writer result but do not propagate it yet: if blastn exited early
        // (e.g. the database was not found) it closes stdin, so the writer sees a benign
        // BrokenPipe. The real cause is the non-zero exit status, so report that first
        // and surface blastn's own stderr message.
        let write_result = writer.join().unwrap();

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            let stderr = stderr.trim();
            if stderr.is_empty() {
                bail!("blastn exited with non-zero status: {}", output.status);
            } else {
                bail!(
                    "blastn exited with non-zero status: {}\n{}",
                    output.status,
                    stderr
                );
            }
        }

        write_result.with_context(|| "Failed to write query to blastn stdin")?;

        // Parse format-6 output from stdout bytes.
        let mut hits_by_index: std::collections::HashMap<usize, Vec<BlastHitFormat6>> =
            std::collections::HashMap::new();
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(output.stdout.as_slice());
        for result in rdr.deserialize::<BlastHitFormat6>() {
            let hit = result?;
            let idx: usize = hit
                .qseqid
                .parse()
                .with_context(|| format!("Malformed query id: {}", hit.qseqid))?;
            hits_by_index.entry(idx).or_default().push(hit);
        }

        let mut results = Vec::with_capacity(baits.len());
        for i in 0..baits.len() {
            let mut hits = hits_by_index.remove(&i).unwrap_or_default();
            // Sort by e-value ascending, then identity descending, then by subject coordinate.
            hits.sort_by(|a, b| {
                a.evalue
                    .partial_cmp(&b.evalue)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| {
                        b.pident
                            .partial_cmp(&a.pident)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .then_with(|| a.sseqid.cmp(&b.sseqid))
                    .then_with(|| a.sstart.cmp(&b.sstart))
            });
            results.push(hits);
        }
        Ok(results)
    }
}

/// Check that the BLAST database exists by looking for index files.
///
/// Checks for `.nsq` / `.nal` (BLAST v4), `.nsi` (BLAST v5), and sharded
/// volume files `.00.nsq` (multi-volume databases such as `nt`).
///
/// When `db_path` is `None`, also reads `$BLASTDB` from the environment and
/// checks there as well.
pub fn database_exists(db: &str, db_path: Option<&Path>) -> bool {
    // Single-volume extensions (v4 and v5).
    let exts = ["nsq", "nal", "nsi"];
    // First-volume extensions for sharded databases.
    let sharded_exts = ["00.nsq", "00.nal"];

    let search_paths: Vec<PathBuf> = {
        let mut paths = Vec::new();
        if let Some(p) = db_path {
            paths.push(p.to_path_buf());
        } else {
            // Fall back to $BLASTDB when no explicit path is given.
            if let Ok(blastdb) = std::env::var("BLASTDB") {
                paths.push(PathBuf::from(blastdb));
            }
            // Also check the current directory / relative path.
            paths.push(PathBuf::new());
        }
        paths
    };

    for dir in &search_paths {
        let in_dir = |ext: &str| -> PathBuf {
            if dir.as_os_str().is_empty() {
                PathBuf::from(format!("{db}.{ext}"))
            } else {
                dir.join(format!("{db}.{ext}"))
            }
        };
        let found = exts.iter().any(|e| in_dir(e).exists())
            || sharded_exts.iter().any(|e| in_dir(e).exists());
        if found {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::skip_if_missing;

    #[test]
    fn test_blast_hit_sort_by_evalue() {
        let mut hits = [
            BlastHitFormat6 {
                qseqid: "0".to_string(),
                sseqid: "chr1".to_string(),
                pident: 100.0,
                length: 100,
                mismatch: 0,
                gapopen: 0,
                qstart: 1,
                qend: 100,
                sstart: 1,
                send: 100,
                evalue: 1e-5,
                bitscore: 200.0,
            },
            BlastHitFormat6 {
                qseqid: "0".to_string(),
                sseqid: "chr2".to_string(),
                pident: 95.0,
                length: 100,
                mismatch: 5,
                gapopen: 0,
                qstart: 1,
                qend: 100,
                sstart: 1,
                send: 100,
                evalue: 1e-2,
                bitscore: 150.0,
            },
        ];
        // already sorted; swap and re-sort
        hits.swap(0, 1);
        hits.sort_by(|a, b| {
            a.evalue
                .partial_cmp(&b.evalue)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        assert_eq!(hits[0].evalue, 1e-5);
        assert_eq!(hits[1].evalue, 1e-2);
    }

    #[test]
    fn test_database_exists_returns_false_for_missing_db() {
        assert!(!database_exists("nonexistent_db_xyz", None));
    }

    #[test]
    fn test_database_exists_returns_false_for_nonexistent_sharded() {
        assert!(!database_exists("nonexistent_sharded_db", None));
    }

    #[test]
    fn test_align_batch_empty_returns_empty_without_blastn() {
        let runner = BlastRunner::new("dummy_db", None, false, 1);
        let results = runner.align_batch(&[]).unwrap();
        assert!(results.is_empty());
    }

    #[test]
    fn test_align_batch_bait_without_sequence_returns_error() {
        let runner = BlastRunner::new("dummy_db", None, false, 1);
        let bait = Bait::new("chr1", 0, 100, "b");
        let result = runner.align_batch(&[bait]);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("no sequence"));
    }

    #[test]
    fn test_database_exists_with_nsq_file() {
        let dir = tempfile::TempDir::new().unwrap();
        std::fs::File::create(dir.path().join("mydb.nsq")).unwrap();
        assert!(database_exists("mydb", Some(dir.path())));
    }

    #[test]
    fn test_database_exists_with_sharded_nsq_file() {
        let dir = tempfile::TempDir::new().unwrap();
        std::fs::File::create(dir.path().join("mydb.00.nsq")).unwrap();
        assert!(database_exists("mydb", Some(dir.path())));
    }

    #[test]
    fn test_database_exists_with_nal_file() {
        let dir = tempfile::TempDir::new().unwrap();
        std::fs::File::create(dir.path().join("mydb.nal")).unwrap();
        assert!(database_exists("mydb", Some(dir.path())));
    }

    #[test]
    fn test_blast_hit_format6_deserialize_from_tsv() {
        let tsv = "query1\tchr1\t100.0\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200.0\n";
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(tsv.as_bytes());
        let hit: BlastHitFormat6 = rdr
            .deserialize::<BlastHitFormat6>()
            .next()
            .unwrap()
            .unwrap();
        assert_eq!(hit.qseqid, "query1");
        assert_eq!(hit.sseqid, "chr1");
        assert!((hit.pident - 100.0).abs() < 1e-9);
        assert_eq!(hit.length, 100);
        assert_eq!(hit.mismatch, 0);
        assert_eq!(hit.gapopen, 0);
        assert_eq!(hit.qstart, 1);
        assert_eq!(hit.qend, 100);
        assert_eq!(hit.sstart, 1);
        assert_eq!(hit.send, 100);
        assert!((hit.evalue - 1e-50).abs() < 1e-60);
        assert!((hit.bitscore - 200.0).abs() < 1e-9);
    }

    #[test]
    fn test_blast_hit_sort_by_pident_when_evalue_equal() {
        // Two hits with equal evalue; the one with higher pident should sort first.
        let mut hits = [
            BlastHitFormat6 {
                qseqid: "0".to_string(),
                sseqid: "chr1".to_string(),
                pident: 90.0,
                length: 100,
                mismatch: 10,
                gapopen: 0,
                qstart: 1,
                qend: 100,
                sstart: 1,
                send: 100,
                evalue: 1e-5,
                bitscore: 150.0,
            },
            BlastHitFormat6 {
                qseqid: "0".to_string(),
                sseqid: "chr2".to_string(),
                pident: 100.0,
                length: 100,
                mismatch: 0,
                gapopen: 0,
                qstart: 1,
                qend: 100,
                sstart: 1,
                send: 100,
                evalue: 1e-5,
                bitscore: 200.0,
            },
        ];
        hits.sort_by(|a, b| {
            a.evalue
                .partial_cmp(&b.evalue)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    b.pident
                        .partial_cmp(&a.pident)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| a.sseqid.cmp(&b.sseqid))
                .then_with(|| a.sstart.cmp(&b.sstart))
        });
        assert!(
            (hits[0].pident - 100.0).abs() < 1e-9,
            "higher pident should sort first"
        );
        assert!((hits[1].pident - 90.0).abs() < 1e-9);
    }

    #[test]
    fn test_blast_hit_sort_by_sseqid_when_evalue_and_pident_equal() {
        let make_hit = |sseqid: &str| BlastHitFormat6 {
            qseqid: "0".to_string(),
            sseqid: sseqid.to_string(),
            pident: 100.0,
            length: 100,
            mismatch: 0,
            gapopen: 0,
            qstart: 1,
            qend: 100,
            sstart: 1,
            send: 100,
            evalue: 1e-5,
            bitscore: 200.0,
        };
        let mut hits = [make_hit("chr3"), make_hit("chr1"), make_hit("chr2")];
        hits.sort_by(|a, b| {
            a.evalue
                .partial_cmp(&b.evalue)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    b.pident
                        .partial_cmp(&a.pident)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| a.sseqid.cmp(&b.sseqid))
                .then_with(|| a.sstart.cmp(&b.sstart))
        });
        assert_eq!(hits[0].sseqid, "chr1");
        assert_eq!(hits[1].sseqid, "chr2");
        assert_eq!(hits[2].sseqid, "chr3");
    }

    #[test]
    fn test_database_exists_with_nsi_file() {
        let dir = tempfile::TempDir::new().unwrap();
        std::fs::File::create(dir.path().join("mydb.nsi")).unwrap();
        assert!(database_exists("mydb", Some(dir.path())));
    }

    #[test]
    fn test_database_exists_with_sharded_nal_file() {
        let dir = tempfile::TempDir::new().unwrap();
        std::fs::File::create(dir.path().join("mydb.00.nal")).unwrap();
        assert!(database_exists("mydb", Some(dir.path())));
    }

    #[test]
    fn test_database_exists_with_blastdb_env() {
        let dir = tempfile::TempDir::new().unwrap();
        std::fs::File::create(dir.path().join("envdb.nsq")).unwrap();
        // Override BLASTDB for this test (safe under --test-threads 1).
        let prev = std::env::var("BLASTDB").ok();
        // SAFETY: test runs with --test-threads 1 (tarpaulin enforces this);
        // no other thread is reading the environment concurrently.
        unsafe {
            std::env::set_var("BLASTDB", dir.path().to_str().unwrap());
        }
        let found = database_exists("envdb", None);
        unsafe {
            match prev {
                Some(v) => std::env::set_var("BLASTDB", v),
                None => std::env::remove_var("BLASTDB"),
            }
        }
        assert!(found);
    }

    #[test]
    fn test_blast_runner_new_stores_fields() {
        let runner = BlastRunner::new("mydb", Some(std::path::PathBuf::from("/tmp")), true, 4);
        assert_eq!(runner.db, "mydb");
        assert_eq!(runner.db_path, Some(std::path::PathBuf::from("/tmp")));
        assert!(runner.dust);
        assert_eq!(runner.num_threads, 4);
    }

    /// Run only when `blastn` and `makeblastdb` are both on $PATH.
    #[test]
    fn test_blast_runner_align_batch_if_available() {
        skip_if_missing!("blastn", "makeblastdb");
        use std::process::Command as StdCommand;
        use tempfile::TempDir;

        let dir = TempDir::new().unwrap();
        let db_prefix = dir.path().join("testdb");

        // 100-bp sequence used as both the database subject and the query.
        let seq: String = "ACGT".repeat(25);
        let seq_len = seq.len() as u64;

        // Write subject FASTA.
        let fasta_path = dir.path().join("subject.fa");
        std::fs::write(&fasta_path, format!(">chr1\n{seq}\n")).unwrap();

        // Build BLAST nucleotide database.
        let mdb = StdCommand::new("makeblastdb")
            .args([
                "-in",
                fasta_path.to_str().unwrap(),
                "-dbtype",
                "nucl",
                "-out",
                db_prefix.to_str().unwrap(),
            ])
            .output()
            .expect("makeblastdb failed to spawn");
        assert!(mdb.status.success(), "makeblastdb exited with failure");

        // Query with the exact same sequence — expect a perfect hit.
        let bait = Bait::with_sequence("chr1", 0, seq_len, "bait0", seq.clone());
        let runner = BlastRunner::new("testdb", Some(dir.path().to_path_buf()), false, 1);

        let results = runner.align_batch(&[bait]).expect("align_batch failed");
        assert_eq!(results.len(), 1);

        let hits = &results[0];
        assert!(
            !hits.is_empty(),
            "expected at least one BLAST hit for perfect match"
        );

        let top = &hits[0];
        assert_eq!(top.sseqid, "chr1");
        assert!((top.pident - 100.0).abs() < 0.01, "expected 100% identity");
        assert_eq!(top.length, seq_len as u32, "expected full-length alignment");
        assert_eq!(top.mismatch, 0);
        assert_eq!(top.gapopen, 0);
        assert!(
            top.evalue < 1e-10,
            "expected tiny e-value, got {}",
            top.evalue
        );
    }
}
