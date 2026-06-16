//! Mappability bedGraph loading and per-bait score extraction.
//!
//! Mappability scores are stored as a bgzf-compressed, tabix-indexed bedGraph
//! where each record assigns a score (0–1) to a genomic interval. Each bait
//! query opens the file, seeks to the relevant block via the `.tbi` index, and
//! expands overlapping records to per-base scores. No genome-wide data is held
//! in memory. Callers apply scores to a metric via
//! [`crate::metrics::apply_mappability`].

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};

use crate::intervals::Bait;
use noodles::core::Region;
use noodles::csi::BinningIndex;
use noodles::tabix;

/// A tabix-indexed bgzf-compressed bedGraph reader for random-access mappability lookup.
///
/// Only a `PathBuf` is stored; each call to [`scores_for_bait`] opens a fresh
/// file handle so the reader is `Send + Sync` and safe for parallel use.
///
/// [`scores_for_bait`]: MappabilityReader::scores_for_bait
#[derive(Debug)]
pub struct MappabilityReader {
    path: PathBuf,
}

impl MappabilityReader {
    /// Open a tabix-indexed bgzf-compressed bedGraph file.
    ///
    /// Expects a `{path}.tbi` index file alongside the compressed file.
    pub fn new(path: &Path) -> Result<Self> {
        if !path.exists() {
            bail!("Mappability file not found: {}", path.display());
        }
        // Build the sibling `.tbi` path from the raw OS bytes, not the lossy
        // `display()` string, so non-UTF-8 filenames resolve to the correct path.
        let tbi = {
            let mut p = path.as_os_str().to_os_string();
            p.push(".tbi");
            PathBuf::from(p)
        };
        if !tbi.exists() {
            bail!(
                "Tabix index not found: {}; run `tabix -p bed {}` to create it",
                tbi.display(),
                path.display()
            );
        }
        Ok(MappabilityReader {
            path: path.to_path_buf(),
        })
    }

    /// Return the set of reference (contig) names present in the tabix index.
    ///
    /// Used to validate bait contigs up front so a missing contig or naming mismatch
    /// (e.g. `chr1` vs `1`) fails fast with a clear message, rather than erroring
    /// mid-run on the first affected bait.
    pub fn reference_names(&self) -> Result<HashSet<String>> {
        let reader = tabix::io::indexed_reader::Builder::default()
            .build_from_path(&self.path)
            .with_context(|| format!("Cannot open tabix file: {}", self.path.display()))?;
        let header = reader
            .index()
            .header()
            .with_context(|| format!("Tabix index has no header: {}", self.path.display()))?;
        Ok(header
            .reference_sequence_names()
            .iter()
            .map(|name| name.to_string())
            .collect())
    }

    /// Extract per-base mappability scores for a given bait interval via tabix query.
    ///
    /// Each bedGraph record overlapping the bait is expanded to per-base scores.
    pub fn scores_for_bait(&self, bait: &Bait) -> Result<Vec<f64>> {
        let region_str = format!("{}:{}-{}", bait.chrom, bait.start + 1, bait.end);
        let region: Region = region_str
            .parse()
            .with_context(|| format!("Invalid region string: {region_str}"))?;

        let mut reader = tabix::io::indexed_reader::Builder::default()
            .build_from_path(&self.path)
            .with_context(|| format!("Cannot open tabix file: {}", self.path.display()))?;

        let query = reader
            .query(&region)
            .with_context(|| format!("Tabix query failed for {region_str}"))?;

        // Initialize all positions to 0.0 which represents no mappability.
        let bait_len = (bait.end - bait.start) as usize;
        let mut scores = vec![0.0_f64; bait_len];
        for result in query {
            let record = result?;
            let line = record.as_ref();
            let fields: Vec<&str> = line.splitn(5, '\t').collect();
            if fields.len() < 4 {
                continue;
            }
            let rec_start: u64 = fields[1]
                .parse()
                .with_context(|| format!("Invalid bedGraph start: {}", fields[1]))?;
            let rec_end: u64 = fields[2]
                .parse()
                .with_context(|| format!("Invalid bedGraph end: {}", fields[2]))?;
            let value_str = fields[3].trim();
            let value: f64 = value_str
                .parse()
                .with_context(|| format!("Invalid bedGraph value: {}", fields[3]))?;

            let overlap_start = rec_start.max(bait.start);
            let overlap_end = rec_end.min(bait.end);
            if overlap_start >= overlap_end {
                continue;
            }
            let idx_start = (overlap_start - bait.start) as usize;
            let idx_end = (overlap_end - bait.start) as usize;
            for s in &mut scores[idx_start..idx_end] {
                *s = value;
            }
        }
        Ok(scores)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_bait(chrom: &str, start: u64, end: u64) -> Bait {
        Bait::new(chrom, start, end, "b")
    }

    /// Compute per-base scores from in-memory bedGraph records (for unit testing).
    ///
    /// Positions not covered by any record are 0.0 (unmappable).
    fn scores_from_records(bait: &Bait, records: &[(&str, u64, u64, f64)]) -> Vec<f64> {
        let bait_len = (bait.end - bait.start) as usize;
        let mut scores = vec![0.0_f64; bait_len];
        for &(chrom, rec_start, rec_end, value) in records {
            if chrom != bait.chrom {
                continue;
            }
            let overlap_start = rec_start.max(bait.start);
            let overlap_end = rec_end.min(bait.end);
            if overlap_start >= overlap_end {
                continue;
            }
            let idx_start = (overlap_start - bait.start) as usize;
            let idx_end = (overlap_end - bait.start) as usize;
            for s in &mut scores[idx_start..idx_end] {
                *s = value;
            }
        }
        scores
    }

    #[test]
    fn test_scores_for_bait_no_overlap_yields_zeros() {
        // No overlapping record → all positions are 0.0 (unmappable).
        let bait = make_bait("chr1", 0, 100);
        let scores = scores_from_records(&bait, &[("chr1", 200, 300, 1.0)]);
        assert_eq!(scores.len(), 100);
        assert!(scores.iter().all(|&s| s == 0.0));
    }

    #[test]
    fn test_scores_for_bait_different_chrom_yields_zeros() {
        let bait = make_bait("chr1", 0, 100);
        let scores = scores_from_records(&bait, &[("chr2", 0, 100, 1.0)]);
        assert_eq!(scores.len(), 100);
        assert!(scores.iter().all(|&s| s == 0.0));
    }

    #[test]
    fn test_scores_for_bait_exact_overlap() {
        // Positions 0,1,2,3 have scores 0.00, 0.25, 0.75, 1.00 (0-based half-open).
        // Bait: chr1 [0, 4) → 4 scores
        let bait = make_bait("chr1", 0, 4);
        let records = [
            ("chr1", 0, 1, 0.00),
            ("chr1", 1, 2, 0.25),
            ("chr1", 2, 3, 0.75),
            ("chr1", 3, 4, 1.00),
        ];
        let scores = scores_from_records(&bait, &records);
        assert_eq!(scores.len(), 4);
        assert_eq!(scores, vec![0.00, 0.25, 0.75, 1.00]);
    }

    #[test]
    fn test_scores_for_bait_partial_overlap() {
        // Bait [5,15): record [0,10) overlaps [5,10) → first 5 bases = 0.5, last 5 = 0.0.
        let bait = make_bait("chr1", 5, 15);
        let scores = scores_from_records(&bait, &[("chr1", 0, 10, 0.5)]);
        assert_eq!(scores.len(), 10);
        assert!(scores[..5].iter().all(|&s| (s - 0.5).abs() < f64::EPSILON));
        assert!(scores[5..].iter().all(|&s| s == 0.0));
    }

    #[test]
    fn test_mappability_reader_new_file_not_found() {
        let result = MappabilityReader::new(std::path::Path::new("/nonexistent/path.bedgraph.gz"));
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("not found"));
    }

    #[test]
    fn test_mappability_reader_new_tbi_not_found() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "placeholder").unwrap();
        let result = MappabilityReader::new(f.path());
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("Tabix index not found")
        );
    }

    /// Path to the committed mappability fixture: `tests/data/mappability.bedgraph.gz` (+ `.tbi`).
    ///
    /// Single record: `chr1  0  200  0.75`.
    fn mappability_fixture() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mappability.bedgraph.gz")
    }

    #[test]
    fn test_mappability_reader_new_success() {
        let reader = MappabilityReader::new(&mappability_fixture());
        assert!(
            reader.is_ok(),
            "MappabilityReader::new failed: {:?}",
            reader.unwrap_err()
        );
    }

    #[test]
    fn test_scores_for_bait_fixture_full_overlap() {
        // Fixture chr1 [0,200) at 0.75; bait [0,10) is fully covered → all 10 scores = 0.75.
        let reader = MappabilityReader::new(&mappability_fixture()).unwrap();
        let bait = make_bait("chr1", 0, 10);
        let scores = reader.scores_for_bait(&bait).unwrap();
        assert_eq!(scores.len(), 10);
        assert!(scores.iter().all(|&s| (s - 0.75).abs() < f64::EPSILON));
    }

    #[test]
    fn test_scores_for_bait_fixture_partial_overlap() {
        // Fixture chr1 [0,200) at 0.75; bait [100,250) extends past the record end.
        // Positions 0..100 of bait (genome [100,200)) = 0.75; positions 100..150 (genome [200,250)) = 0.0.
        let reader = MappabilityReader::new(&mappability_fixture()).unwrap();
        let bait = make_bait("chr1", 100, 250);
        let scores = reader.scores_for_bait(&bait).unwrap();
        assert_eq!(scores.len(), 150);
        assert!(
            scores[..100]
                .iter()
                .all(|&s| (s - 0.75).abs() < f64::EPSILON)
        );
        assert!(scores[100..].iter().all(|&s| s == 0.0));
    }

    #[test]
    fn test_scores_for_bait_fixture_no_overlap_returns_zeros() {
        // Fixture chr1 [0,200); bait [300,310) has no overlapping records → all zeros.
        let reader = MappabilityReader::new(&mappability_fixture()).unwrap();
        let bait = make_bait("chr1", 300, 310);
        let scores = reader.scores_for_bait(&bait).unwrap();
        assert!(scores.iter().all(|&s| s == 0.0));
    }
}
