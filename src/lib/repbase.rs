//! RepBase repeat-feature annotation via tabix-indexed BED overlap.
//!
//! A bgzf-compressed, tabix-indexed BED file of RepBase features is queried
//! per bait. For each bait, overlapping feature names are collected, sorted,
//! and returned as a `Vec<String>` for the caller to apply via
//! [`crate::metrics::apply_repbase`]. No genome-wide data is held in memory.

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};

use crate::intervals::Bait;
use noodles::core::Region;
use noodles::csi::BinningIndex;
use noodles::tabix;

/// A tabix-indexed bgzf-compressed BED reader for random-access RepBase lookup.
///
/// Only a `PathBuf` is stored; each call to [`overlapping_features`] opens a
/// fresh file handle so the reader is `Send + Sync` and safe for parallel use.
///
/// [`overlapping_features`]: RepBaseReader::overlapping_features
#[derive(Debug)]
pub struct RepBaseReader {
    path: PathBuf,
}

impl RepBaseReader {
    /// Open a tabix-indexed bgzf-compressed BED file.
    ///
    /// Expects a `{path}.tbi` index file alongside the compressed file.
    pub fn new(path: &Path) -> Result<Self> {
        if !path.exists() {
            bail!("RepBase file not found: {}", path.display());
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
        Ok(RepBaseReader {
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

    /// Return a sorted list of feature names overlapping `bait`.
    pub fn overlapping_features(&self, bait: &Bait) -> Result<Vec<String>> {
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

        let mut names = Vec::new();
        for result in query {
            let record = result?;
            let line = record.as_ref();
            let fields: Vec<&str> = line.splitn(5, '\t').collect();
            if fields.len() < 4 {
                continue;
            }
            names.push(fields[3].trim().to_string());
        }
        names.sort();
        names.dedup();
        Ok(names)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_bait(chrom: &str, start: u64, end: u64) -> Bait {
        Bait::new(chrom, start, end, "b")
    }

    /// Return sorted feature names overlapping `bait` from in-memory records (for unit testing).
    fn features_from_records(bait: &Bait, records: &[(&str, u64, u64, &str)]) -> Vec<String> {
        let mut names = Vec::new();
        for &(chrom, rec_start, rec_end, name) in records {
            if chrom != bait.chrom {
                continue;
            }
            // Tabix-style overlap: any overlap between [rec_start, rec_end) and [bait.start, bait.end)
            if rec_start < bait.end && rec_end > bait.start {
                names.push(name.to_string());
            }
        }
        names.sort();
        names.dedup();
        names
    }

    #[test]
    fn test_no_overlap_returns_empty() {
        // TAR1 at [0,1) is 1bp left of bait [1,6) → no overlap
        let bait = make_bait("chr1", 1, 6);
        let features = features_from_records(&bait, &[("chr1", 0, 1, "TAR1")]);
        assert!(features.is_empty());
    }

    #[test]
    fn test_overlap_returns_name() {
        let bait = make_bait("chr1", 1, 6);
        let features = features_from_records(&bait, &[("chr1", 1, 6, "(TAACCC)n")]);
        assert_eq!(features, vec!["(TAACCC)n"]);
    }

    #[test]
    fn test_repbase_overlap_case() {
        // Bait [1,6), features:
        //   TAR1 [0,1) → no overlap (1 bp to the left)
        //   (TAACCC)n [1,6) → overlap
        //   L1MC5a [4,5) → overlap (within bait)
        //   MIR3 [6,7) → no overlap (1 bp to the right)
        let bait = make_bait("chr1", 1, 6);
        let records = [
            ("chr1", 0, 1, "TAR1"),
            ("chr1", 1, 6, "(TAACCC)n"),
            ("chr1", 4, 5, "L1MC5a"),
            ("chr1", 6, 7, "MIR3"),
        ];
        let features = features_from_records(&bait, &records);
        assert_eq!(features, vec!["(TAACCC)n", "L1MC5a"]);
    }

    #[test]
    fn test_features_sorted() {
        let bait = make_bait("chr1", 0, 100);
        let records = [
            ("chr1", 0, 100, "Zzz"),
            ("chr1", 0, 100, "Aaa"),
            ("chr1", 0, 100, "Mmm"),
        ];
        let features = features_from_records(&bait, &records);
        assert_eq!(features, vec!["Aaa", "Mmm", "Zzz"]);
    }

    #[test]
    fn test_wrong_chrom_returns_empty() {
        let bait = make_bait("chr2", 0, 100);
        let features = features_from_records(&bait, &[("chr1", 0, 100, "Alu")]);
        assert!(features.is_empty());
    }

    #[test]
    fn test_features_deduped() {
        let bait = make_bait("chr1", 0, 100);
        // Two records with the same name → should appear only once.
        let records = [("chr1", 0, 50, "Alu"), ("chr1", 50, 100, "Alu")];
        let features = features_from_records(&bait, &records);
        assert_eq!(features, vec!["Alu"]);
    }

    #[test]
    fn test_repbase_reader_new_file_not_found() {
        let result = RepBaseReader::new(std::path::Path::new("/nonexistent/path.bed.gz"));
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("not found"));
    }

    #[test]
    fn test_repbase_reader_new_tbi_not_found() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "placeholder").unwrap();
        let result = RepBaseReader::new(f.path());
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("Tabix index not found")
        );
    }

    #[test]
    fn test_repbase_reader_new_success() {
        let reader = RepBaseReader::new(&repbase_fixture());
        assert!(
            reader.is_ok(),
            "RepBaseReader::new failed: {:?}",
            reader.unwrap_err()
        );
    }

    /// Path to the committed repbase fixture: `tests/data/repbase.bed.gz` (+ `.tbi`).
    ///
    /// Records:
    /// ```text
    /// chr1  0    50   Alu
    /// chr1  0    100  L1
    /// chr1  50   100  Alu   ← duplicate of the first Alu
    /// chr1  500  600  MIR3
    /// ```
    fn repbase_fixture() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/repbase.bed.gz")
    }

    #[test]
    fn test_overlapping_features_fixture_overlapping_sorted_deduped() {
        // Bait [0,100) overlaps Alu (twice) and L1 → sorted+deduped: ["Alu", "L1"]
        let reader = RepBaseReader::new(&repbase_fixture()).unwrap();
        let bait = make_bait("chr1", 0, 100);
        let features = reader.overlapping_features(&bait).unwrap();
        assert_eq!(features, vec!["Alu", "L1"]);
    }

    #[test]
    fn test_overlapping_features_fixture_single_hit() {
        // Bait [500,600) overlaps only MIR3
        let reader = RepBaseReader::new(&repbase_fixture()).unwrap();
        let bait = make_bait("chr1", 500, 600);
        let features = reader.overlapping_features(&bait).unwrap();
        assert_eq!(features, vec!["MIR3"]);
    }

    #[test]
    fn test_overlapping_features_fixture_no_overlap() {
        // Bait [200,300) has no overlapping records
        let reader = RepBaseReader::new(&repbase_fixture()).unwrap();
        let bait = make_bait("chr1", 200, 300);
        let features = reader.overlapping_features(&bait).unwrap();
        assert!(features.is_empty());
    }

    #[test]
    fn test_overlapping_features_fixture_partial_overlap_left() {
        // Bait [40,60) clips both Alu records and L1 → sorted+deduped: ["Alu", "L1"]
        let reader = RepBaseReader::new(&repbase_fixture()).unwrap();
        let bait = make_bait("chr1", 40, 60);
        let features = reader.overlapping_features(&bait).unwrap();
        assert_eq!(features, vec!["Alu", "L1"]);
    }
}
