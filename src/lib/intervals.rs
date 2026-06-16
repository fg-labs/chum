//! BED, Interval List, and FASTA bait parsing with auto-format detection.
//!
//! Format detection works by peeking at the first byte of the input source
//! without consuming it ([`BufRead::fill_buf`]):
//!   - `@` → Interval List (Picard SAM-header style)
//!   - `>` → FASTA (sequence pre-loaded, coordinates extracted from header)
//!   - Anything else → BED
//!
//! All parsing functions come in two flavors: a path-based convenience
//! wrapper and a `_reader` variant that accepts any [`BufRead`] source,
//! enabling streaming use from stdin, compressed streams, or test cursors.
//!
//! BED records use 0-based half-open coordinates. Interval List records use
//! 1-based closed coordinates and are converted to 0-based on parse. FASTA
//! headers are expected to contain a `chrom:start-end` token in UCSC style
//! (1-based fully closed); it is converted to 0-based half-open on parse.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};
use bio::io::fasta;
use noodles::sam;

/// A genomic interval representing a single bait.
#[derive(Debug, Clone, PartialEq)]
pub struct Bait {
    /// Reference sequence name (chromosome).
    pub chrom: String,
    /// 0-based start position (inclusive).
    pub start: u64,
    /// 0-based end position (exclusive).
    pub end: u64,
    /// Bait name from the interval name field.
    pub name: String,
    /// Strand (`+`, `-`, or `.`), when available from the input format.
    pub strand: Option<char>,
    /// Nucleotide sequence, populated after reference FASTA lookup.
    pub sequence: Option<String>,
}

impl Bait {
    /// Create a new bait from 0-based half-open coordinates.
    pub fn new(chrom: impl Into<String>, start: u64, end: u64, name: impl Into<String>) -> Self {
        Bait {
            chrom: chrom.into(),
            start,
            end,
            name: name.into(),
            strand: None,
            sequence: None,
        }
    }

    /// Create a bait with a pre-loaded sequence.
    pub fn with_sequence(
        chrom: impl Into<String>,
        start: u64,
        end: u64,
        name: impl Into<String>,
        sequence: impl Into<String>,
    ) -> Self {
        Bait {
            chrom: chrom.into(),
            start,
            end,
            name: name.into(),
            strand: None,
            sequence: Some(sequence.into()),
        }
    }

    /// Length of the bait in bases.
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// `true` when the bait has zero length.
    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }

    /// `true` when the bait has real genomic coordinates.
    ///
    /// Returns `false` for FASTA-input baits whose headers contained no
    /// parseable `chrom:start-end` token; those baits receive the synthetic
    /// coordinates `unknown:0-<seq_len>` and cannot be used for tabix or
    /// target-centering queries.
    pub fn is_located(&self) -> bool {
        self.chrom != "unknown"
    }
}

/// Supported interval file formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntervalFormat {
    /// BED format (0-based half-open, tab-separated).
    Bed,
    /// Picard/GATK Interval List format (SAM header + 1-based closed records).
    IntervalList,
    /// FASTA format — coordinates extracted from `>name chrom:start-end` headers.
    Fasta,
}

/// Detect the interval format by peeking at the first byte without consuming it.
///
/// Uses [`BufRead::fill_buf`] so the reader stays positioned at byte 0.
///
/// - `@` → Interval List
/// - `>` → FASTA
/// - otherwise → BED
pub fn detect_format_reader<R: BufRead>(reader: &mut R) -> Result<IntervalFormat> {
    let buf = reader
        .fill_buf()
        .with_context(|| "I/O error peeking interval source")?;
    if buf.is_empty() {
        bail!("Interval source is empty");
    }
    Ok(match buf[0] {
        b'@' => IntervalFormat::IntervalList,
        b'>' => IntervalFormat::Fasta,
        _ => IntervalFormat::Bed,
    })
}

/// Detect the interval format by peeking at the first byte of the file.
///
/// - `@` → Interval List
/// - `>` → FASTA
/// - otherwise → BED
pub fn detect_format(path: &Path) -> Result<IntervalFormat> {
    let file =
        File::open(path).with_context(|| format!("Cannot open baits file: {}", path.display()))?;
    detect_format_reader(&mut BufReader::new(file))
        .with_context(|| format!("in {}", path.display()))
}

/// Parse baits from a [`BufRead`] source, auto-detecting the format.
pub fn parse_intervals_reader<R: BufRead>(mut reader: R) -> Result<Vec<Bait>> {
    let format = detect_format_reader(&mut reader)?;
    match format {
        IntervalFormat::Bed => parse_bed_reader(reader),
        IntervalFormat::IntervalList => parse_interval_list_reader(reader),
        IntervalFormat::Fasta => parse_fasta_baits_reader(reader),
    }
}

/// Parse baits from a BED, Interval List, or FASTA file, auto-detecting the format.
pub fn parse_intervals(path: &Path) -> Result<Vec<Bait>> {
    let file =
        File::open(path).with_context(|| format!("Cannot open baits file: {}", path.display()))?;
    parse_intervals_reader(BufReader::new(file))
}

/// Parse baits from a BED file (0-based, half-open).
///
/// Lines beginning with `#` or `track` or `browser` are skipped. Only the first
/// four columns (chrom, start, end, name) are required; additional columns are ignored.
/// If the name column is absent the interval coordinates are used as the name.
pub fn parse_bed(path: &Path) -> Result<Vec<Bait>> {
    let file =
        File::open(path).with_context(|| format!("Cannot open BED file: {}", path.display()))?;
    parse_bed_reader(BufReader::new(file))
}

/// Parse BED records from any `BufRead` source.
pub fn parse_bed_reader<R: BufRead>(reader: R) -> Result<Vec<Bait>> {
    let mut baits = Vec::new();
    for (lineno, line_result) in reader.lines().enumerate() {
        let line = line_result.with_context(|| format!("I/O error on BED line {}", lineno + 1))?;
        let line = line.trim();
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let fields: Vec<&str> = line.splitn(7, '\t').collect();
        if fields.len() < 3 {
            bail!(
                "BED line {} has fewer than 3 fields: {:?}",
                lineno + 1,
                line
            );
        }
        let chrom = fields[0].to_string();
        let start: u64 = fields[1]
            .parse()
            .with_context(|| format!("BED line {}: invalid start '{}'", lineno + 1, fields[1]))?;
        let end: u64 = fields[2]
            .parse()
            .with_context(|| format!("BED line {}: invalid end '{}'", lineno + 1, fields[2]))?;
        if start > end {
            bail!("BED line {}: start ({}) > end ({})", lineno + 1, start, end);
        }
        let name = if fields.len() >= 4 && !fields[3].is_empty() {
            fields[3].to_string()
        } else {
            format!("{}:{}-{}", chrom, start, end)
        };
        let strand = fields
            .get(5)
            .and_then(|s| s.chars().next())
            .filter(|c| matches!(c, '+' | '-' | '.'));
        let mut bait = Bait::new(chrom, start, end, name);
        bait.strand = strand;
        baits.push(bait);
    }
    Ok(baits)
}

/// Parse baits from a Picard/GATK Interval List file (1-based closed coordinates).
///
/// SAM header lines (`@`-prefixed) are consumed with [`noodles::sam`]; data lines have
/// the format `chrom\tstart\tend\tstrand\tname` and are converted to 0-based half-open.
pub fn parse_interval_list(path: &Path) -> Result<Vec<Bait>> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open Interval List file: {}", path.display()))?;
    parse_interval_list_reader(BufReader::new(file))
}

/// Parse Interval List records from any `BufRead` source.
///
/// Uses [`noodles::sam`] to consume all `@`-prefixed SAM header lines, leaving
/// the reader positioned at the first data line.
pub fn parse_interval_list_reader<R: BufRead>(mut reader: R) -> Result<Vec<Bait>> {
    // Consume the SAM/Picard header (@-prefixed lines); noodles peeks via fill_buf
    // and stops without consuming the first non-@ byte, so reader stays positioned
    // at the start of the first data line.
    {
        let mut sam_reader = sam::io::Reader::new(&mut reader);
        sam_reader
            .read_header()
            .with_context(|| "Error reading Interval List SAM header")?;
    }
    let mut baits = Vec::new();
    for (lineno, line_result) in reader.lines().enumerate() {
        let line = line_result
            .with_context(|| format!("I/O error on Interval List line {}", lineno + 1))?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.splitn(6, '\t').collect();
        if fields.len() < 5 {
            bail!(
                "Interval List line {} has fewer than 5 fields: {:?}",
                lineno + 1,
                line
            );
        }
        let chrom = fields[0].to_string();
        // Interval List is 1-based closed; convert to 0-based half-open.
        let start_1: u64 = fields[1].parse().with_context(|| {
            format!(
                "Interval List line {}: invalid start '{}'",
                lineno + 1,
                fields[1]
            )
        })?;
        let end_1: u64 = fields[2].parse().with_context(|| {
            format!(
                "Interval List line {}: invalid end '{}'",
                lineno + 1,
                fields[2]
            )
        })?;
        if start_1 == 0 {
            bail!(
                "Interval List line {}: start position is 0; Interval List coordinates are 1-based",
                lineno + 1
            );
        }
        if start_1 > end_1 {
            bail!(
                "Interval List line {}: start ({}) > end ({})",
                lineno + 1,
                start_1,
                end_1
            );
        }
        let start = start_1 - 1;
        let end = end_1;
        let strand = fields[3]
            .chars()
            .next()
            .filter(|c| matches!(c, '+' | '-' | '.'));
        let name = fields[4].to_string();
        let mut bait = Bait::new(chrom, start, end, name);
        bait.strand = strand;
        baits.push(bait);
    }
    Ok(baits)
}

/// Parse baits from a FASTA file.
///
/// Each record's sequence is stored in `bait.sequence`. Genomic coordinates
/// are extracted from the header by scanning whitespace- or pipe-separated
/// tokens for a `chrom:start-end` pattern in UCSC style (1-based fully
/// closed), which is converted to 0-based half-open on parse.
/// If no coordinate token is found, `chrom` is set to `"unknown"`, `start` to
/// `0`, and `end` to the sequence length, so metrics that don't require
/// coordinates (e.g. GC content, homopolymers, RNAFold) still work.
pub fn parse_fasta_baits(path: &Path) -> Result<Vec<Bait>> {
    let file =
        File::open(path).with_context(|| format!("Cannot open FASTA file: {}", path.display()))?;
    parse_fasta_baits_reader(BufReader::new(file))
}

/// Parse FASTA baits from any `BufRead` source.
pub fn parse_fasta_baits_reader<R: BufRead>(reader: R) -> Result<Vec<Bait>> {
    let fasta_reader = fasta::Reader::new(reader);
    let mut baits = Vec::new();
    for result in fasta_reader.records() {
        let record = result.with_context(|| "Error reading FASTA record")?;
        let id = record.id();
        let desc = record.desc().unwrap_or("");
        let seq = String::from_utf8_lossy(record.seq()).into_owned();
        let seq_len = seq.len() as u64;

        // Collect tokens from id and description, splitting on whitespace and '|'.
        let coord = id
            .split('|')
            .chain(desc.split(|c: char| c.is_whitespace() || c == '|'))
            .find_map(parse_coord_token);

        let (chrom, start, end, name) = match coord {
            Some((chrom, start, end)) => (chrom, start, end, id.to_string()),
            None => ("unknown".to_string(), 0u64, seq_len, id.to_string()),
        };

        baits.push(Bait::with_sequence(chrom, start, end, name, seq));
    }
    Ok(baits)
}

/// Try to parse a `chrom:start-end` coordinate from a single token.
///
/// The coordinate is interpreted as UCSC style (1-based fully closed) and
/// converted to 0-based half-open on return. Handles optional trailing
/// characters after the end coordinate (e.g. `chr1:100-220_extra`).
/// Returns `(chrom, start, end)` in 0-based half-open, or `None`.
fn parse_coord_token(token: &str) -> Option<(String, u64, u64)> {
    let colon = token.find(':')?;
    let chrom = token[..colon].to_string();
    if chrom.is_empty() {
        return None;
    }
    let rest = &token[colon + 1..];
    let dash = rest.find('-')?;
    let start_str = &rest[..dash];
    // End may be followed by non-digit characters; take only the digit prefix.
    let end_str: String = rest[dash + 1..]
        .chars()
        .take_while(|c| c.is_ascii_digit())
        .collect();
    // Parse as 1-based closed (UCSC); convert to 0-based half-open.
    let start_1based = start_str.parse::<u64>().ok()?;
    let end = end_str.parse::<u64>().ok()?;
    if start_1based == 0 || start_1based > end {
        return None;
    }
    Some((chrom, start_1based - 1, end))
}

/// A target interval for grouping bait metrics.
#[derive(Debug, Clone, PartialEq)]
pub struct Target {
    /// Reference sequence name.
    pub chrom: String,
    /// 0-based start position (inclusive).
    pub start: u64,
    /// 0-based end position (exclusive).
    pub end: u64,
    /// Target name.
    pub name: String,
}

impl Target {
    /// Create a target from 0-based half-open coordinates.
    pub fn new(chrom: impl Into<String>, start: u64, end: u64, name: impl Into<String>) -> Self {
        Target {
            chrom: chrom.into(),
            start,
            end,
            name: name.into(),
        }
    }

    /// Length of the target in bases.
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    /// `true` when the target has zero length.
    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }

    /// Return a padded copy of this target (clamped to 0 on the left).
    pub fn padded(&self, padding: u32) -> Target {
        Target {
            chrom: self.chrom.clone(),
            start: self.start.saturating_sub(padding as u64),
            end: self.end + padding as u64,
            name: self.name.clone(),
        }
    }
}

/// Parse targets from a [`BufRead`] source, auto-detecting the format.
///
/// FASTA is not a valid target format — targets must have explicit genomic
/// coordinates in BED or Interval List form.
pub fn parse_targets_reader<R: BufRead>(mut reader: R) -> Result<Vec<Target>> {
    let format = detect_format_reader(&mut reader)?;
    let baits = match format {
        IntervalFormat::Bed => parse_bed_reader(reader)?,
        IntervalFormat::IntervalList => parse_interval_list_reader(reader)?,
        IntervalFormat::Fasta => {
            bail!("Target source appears to be FASTA; targets must be BED or Interval List")
        }
    };
    Ok(baits
        .into_iter()
        .map(|b| Target::new(b.chrom, b.start, b.end, b.name))
        .collect())
}

/// Parse targets from a BED or Interval List file, auto-detecting the format.
///
/// FASTA is not a valid target format — targets must have explicit genomic
/// coordinates in BED or Interval List form.
pub fn parse_targets(path: &Path) -> Result<Vec<Target>> {
    let file =
        File::open(path).with_context(|| format!("Cannot open target file: {}", path.display()))?;
    parse_targets_reader(BufReader::new(file)).with_context(|| format!("in {}", path.display()))
}

/// Compute the bait-centering metric for a bait relative to a target.
///
/// Measures how centered the bait–target overlap is within the **bait**:
/// 1.0 means the overlap is perfectly centred in the bait; the value falls
/// toward 0.0 as the overlap shifts to one edge of the bait.
///
/// This equation handles both ends of the size spectrum correctly:
/// - **Small targets (SNPs):** the overlap ≈ the target itself, so the score
///   reflects how well the probe centers on the variant site.
/// - **Large targets (wider than the bait):** the overlap equals the full bait,
///   so every bait that is completely covered by the target scores 1.0
///   regardless of where it sits within the target.
///
/// Returns 0.0 when there is no overlap.
pub fn target_centering(bait: &Bait, target: &Target) -> f64 {
    if bait.chrom != target.chrom {
        return 0.0;
    }
    let overlap_start = bait.start.max(target.start) as f64;
    let overlap_end = bait.end.min(target.end) as f64;
    if overlap_start >= overlap_end {
        return 0.0;
    }
    let overlap_center = (overlap_start + overlap_end) / 2.0;
    let bait_center = (bait.start + bait.end) as f64 / 2.0;
    let bait_len = (bait.end - bait.start) as f64;
    1.0 - (2.0 * (overlap_center - bait_center).abs()) / bait_len
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_bait_len() {
        let b = Bait::new("chr1", 100, 200, "b1");
        assert_eq!(b.len(), 100);
    }

    #[test]
    fn test_bait_is_empty() {
        let b = Bait::new("chr1", 100, 100, "b1");
        assert!(b.is_empty());
    }

    #[test]
    fn test_detect_format_reader_bed() {
        let mut cursor = Cursor::new("chr1\t0\t100\tbait1\n");
        assert_eq!(
            detect_format_reader(&mut cursor).unwrap(),
            IntervalFormat::Bed
        );
        // Reader must still be at the start after peeking.
        let mut line = String::new();
        cursor.read_line(&mut line).unwrap();
        assert!(line.starts_with("chr1"));
    }

    #[test]
    fn test_detect_format_reader_interval_list() {
        let mut cursor = Cursor::new("@HD\tVN:1.6\n");
        assert_eq!(
            detect_format_reader(&mut cursor).unwrap(),
            IntervalFormat::IntervalList
        );
    }

    #[test]
    fn test_detect_format_reader_fasta() {
        let mut cursor = Cursor::new(">seq1\nACGT\n");
        assert_eq!(
            detect_format_reader(&mut cursor).unwrap(),
            IntervalFormat::Fasta
        );
    }

    #[test]
    fn test_detect_format_reader_empty_is_error() {
        let mut cursor = Cursor::new("");
        assert!(detect_format_reader(&mut cursor).is_err());
    }

    #[test]
    fn test_detect_format_bed_from_tmpfile() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t0\t100\tbait1").unwrap();
        assert_eq!(detect_format(f.path()).unwrap(), IntervalFormat::Bed);
    }

    #[test]
    fn test_detect_format_interval_list_from_tmpfile() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "@HD\tVN:1.6").unwrap();
        writeln!(f, "chr1\t1\t100\t+\tbait1").unwrap();
        assert_eq!(
            detect_format(f.path()).unwrap(),
            IntervalFormat::IntervalList
        );
    }

    #[test]
    fn test_parse_intervals_reader_bed() {
        let input = "chr1\t100\t200\tbait1\n";
        let baits = parse_intervals_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits.len(), 1);
        assert_eq!(baits[0], Bait::new("chr1", 100, 200, "bait1"));
    }

    #[test]
    fn test_parse_intervals_reader_interval_list() {
        let input = "@HD\tVN:1.6\nchr1\t101\t200\t+\tbait1\n";
        let baits = parse_intervals_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits.len(), 1);
        assert_eq!(baits[0].chrom, "chr1");
        assert_eq!(baits[0].start, 100);
        assert_eq!(baits[0].end, 200);
    }

    #[test]
    fn test_parse_intervals_reader_fasta() {
        let input = ">bait1 chr1:100-200\nACGT\n";
        let baits = parse_intervals_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits.len(), 1);
        assert_eq!(baits[0].chrom, "chr1");
    }

    #[test]
    fn test_parse_bed_basic() {
        let input = "chr1\t100\t200\tbait_1\nchr2\t300\t400\tbait_2\n";
        let baits = parse_bed_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits.len(), 2);
        assert_eq!(baits[0], Bait::new("chr1", 100, 200, "bait_1"));
        assert_eq!(baits[1], Bait::new("chr2", 300, 400, "bait_2"));
    }

    #[test]
    fn test_parse_bed_skips_comment_lines() {
        let input = "# comment\nchr1\t0\t10\tb1\n";
        let baits = parse_bed_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits.len(), 1);
    }

    #[test]
    fn test_parse_bed_skips_track_lines() {
        let input = "track name=foo\nchr1\t0\t10\tb1\n";
        let baits = parse_bed_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits.len(), 1);
    }

    #[test]
    fn test_parse_bed_skips_browser_lines() {
        let input = "browser position chr1:1-100\nchr1\t0\t10\tb1\n";
        let baits = parse_bed_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits.len(), 1);
    }

    #[test]
    fn test_parse_bed_three_col_uses_coordinates_as_name() {
        let input = "chr1\t100\t200\n";
        let baits = parse_bed_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits[0].name, "chr1:100-200");
    }

    #[test]
    fn test_parse_bed_extra_columns_ignored() {
        let input = "chr1\t100\t200\tbait\t0\t+\n";
        let baits = parse_bed_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits.len(), 1);
        assert_eq!(baits[0].name, "bait");
    }

    #[test]
    fn test_parse_bed_empty_input() {
        let baits = parse_bed_reader(Cursor::new("")).unwrap();
        assert!(baits.is_empty());
    }

    #[test]
    fn test_parse_bed_error_on_too_few_fields() {
        let result = parse_bed_reader(Cursor::new("chr1\t100\n"));
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_interval_list_basic() {
        let input = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nchr1\t101\t200\t+\tbait_1\n";
        let baits = parse_interval_list_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits.len(), 1);
        // 1-based [101, 200] → 0-based [100, 200)
        let mut expected = Bait::new("chr1", 100, 200, "bait_1");
        expected.strand = Some('+');
        assert_eq!(baits[0], expected);
    }

    #[test]
    fn test_parse_interval_list_converts_to_zero_based() {
        // Interval List: 1-based closed [1, 10] → 0-based half-open [0, 10)
        let input = "@HD\tVN:1.6\nchr1\t1\t10\t+\tbait\n";
        let baits = parse_interval_list_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits[0].start, 0);
        assert_eq!(baits[0].end, 10);
        assert_eq!(baits[0].len(), 10);
    }

    #[test]
    fn test_parse_interval_list_skips_at_lines() {
        let input = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nchr1\t1\t10\t+\tb\n";
        let baits = parse_interval_list_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits.len(), 1);
    }

    #[test]
    fn test_parse_interval_list_no_header() {
        // Interval List without any @ header lines should still parse.
        let input = "chr1\t1\t10\t+\tbait\n";
        let baits = parse_interval_list_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits.len(), 1);
        assert_eq!(baits[0].start, 0);
        assert_eq!(baits[0].end, 10);
    }

    #[test]
    fn test_parse_fasta_baits_reader_with_coords() {
        // UCSC 1-based closed [100, 220] → 0-based half-open [99, 220)
        let input = ">bait1 chr1:100-220\nACGTACGT\n";
        let baits = parse_fasta_baits_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits.len(), 1);
        assert_eq!(baits[0].chrom, "chr1");
        assert_eq!(baits[0].start, 99);
        assert_eq!(baits[0].end, 220);
        assert_eq!(baits[0].sequence, Some("ACGTACGT".to_string()));
    }

    #[test]
    fn test_parse_fasta_baits_reader_no_coords() {
        let input = ">just_a_name\nAAAA\n";
        let baits = parse_fasta_baits_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits[0].chrom, "unknown");
        assert_eq!(baits[0].start, 0);
        assert_eq!(baits[0].end, 4);
    }

    #[test]
    fn test_parse_coord_token_basic() {
        // UCSC 1-based closed [100, 220] → 0-based half-open [99, 220)
        assert_eq!(
            parse_coord_token("chr1:100-220"),
            Some(("chr1".to_string(), 99, 220))
        );
    }

    #[test]
    fn test_parse_coord_token_with_trailing_text() {
        assert_eq!(
            parse_coord_token("chr1:100-220_extra"),
            Some(("chr1".to_string(), 99, 220))
        );
    }

    #[test]
    fn test_parse_coord_token_no_colon_returns_none() {
        assert!(parse_coord_token("chr1_100_220").is_none());
    }

    #[test]
    fn test_parse_coord_token_invalid_start_returns_none() {
        assert!(parse_coord_token("chr1:abc-220").is_none());
    }

    #[test]
    fn test_parse_fasta_baits_with_coords_in_id() {
        // UCSC 1-based closed [100, 220] → 0-based half-open [99, 220)
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, ">chr1:100-220").unwrap();
        writeln!(f, "ACGTACGTACGT").unwrap();
        let baits = parse_fasta_baits(f.path()).unwrap();
        assert_eq!(baits.len(), 1);
        assert_eq!(baits[0].chrom, "chr1");
        assert_eq!(baits[0].start, 99);
        assert_eq!(baits[0].end, 220);
        assert_eq!(baits[0].name, "chr1:100-220");
        assert_eq!(baits[0].sequence, Some("ACGTACGTACGT".to_string()));
    }

    #[test]
    fn test_parse_fasta_baits_coords_in_description() {
        // UCSC 1-based closed [100, 220] → 0-based half-open [99, 220)
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, ">bait_name chr1:100-220").unwrap();
        writeln!(f, "AAAA").unwrap();
        let baits = parse_fasta_baits(f.path()).unwrap();
        assert_eq!(baits[0].chrom, "chr1");
        assert_eq!(baits[0].start, 99);
        assert_eq!(baits[0].end, 220);
        assert_eq!(baits[0].name, "bait_name");
    }

    #[test]
    fn test_parse_fasta_baits_no_coords_uses_synthetic() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, ">just_a_name").unwrap();
        writeln!(f, "ACGTACGT").unwrap();
        let baits = parse_fasta_baits(f.path()).unwrap();
        assert_eq!(baits[0].chrom, "unknown");
        assert_eq!(baits[0].start, 0);
        assert_eq!(baits[0].end, 8);
        assert_eq!(baits[0].name, "just_a_name");
    }

    #[test]
    fn test_parse_fasta_baits_pipe_separated_coord() {
        // UCSC 1-based closed [500, 620] → 0-based half-open [499, 620)
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, ">bait_name|chr1:500-620").unwrap();
        writeln!(f, "TTTTCCCC").unwrap();
        let baits = parse_fasta_baits(f.path()).unwrap();
        assert_eq!(baits[0].chrom, "chr1");
        assert_eq!(baits[0].start, 499);
        assert_eq!(baits[0].end, 620);
    }

    #[test]
    fn test_parse_targets_reader_bed() {
        let input = "chr1\t100\t200\ttarget1\n";
        let targets = parse_targets_reader(Cursor::new(input)).unwrap();
        assert_eq!(targets.len(), 1);
        assert_eq!(targets[0], Target::new("chr1", 100, 200, "target1"));
    }

    #[test]
    fn test_parse_targets_reader_interval_list() {
        let input = "@HD\tVN:1.6\nchr1\t101\t200\t+\ttarget1\n";
        let targets = parse_targets_reader(Cursor::new(input)).unwrap();
        assert_eq!(targets.len(), 1);
        assert_eq!(targets[0], Target::new("chr1", 100, 200, "target1"));
    }

    #[test]
    fn test_parse_targets_reader_fasta_is_error() {
        let input = ">seq\nACGT\n";
        assert!(parse_targets_reader(Cursor::new(input)).is_err());
    }

    #[test]
    fn test_target_centering_small_target_centered() {
        // SNP-like target smaller than bait, centred → 1.0
        let bait = Bait::new("chr1", 100, 220, "b"); // center = 160, len = 120
        let target = Target::new("chr1", 155, 165, "snp"); // overlap_center = 160
        assert!((target_centering(&bait, &target) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_target_centering_large_target_bait_fully_covered() {
        // Target larger than bait — bait sits anywhere inside it → always 1.0
        let target = Target::new("chr1", 0, 600, "big");

        let bait_left = Bait::new("chr1", 0, 120, "left");
        assert!((target_centering(&bait_left, &target) - 1.0).abs() < 1e-9);

        let bait_right = Bait::new("chr1", 480, 600, "right");
        assert!((target_centering(&bait_right, &target) - 1.0).abs() < 1e-9);

        let bait_mid = Bait::new("chr1", 240, 360, "mid");
        assert!((target_centering(&bait_mid, &target) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_target_centering_partial_overlap_left_edge() {
        // Target clips only the left 20 bp of a 120 bp bait.
        // overlap [100,120), overlap_center=110, bait_center=160, bait_len=120
        // 1 - 2*50/120 ≈ 0.167
        let bait = Bait::new("chr1", 100, 220, "b");
        let target = Target::new("chr1", 80, 120, "t");
        let c = target_centering(&bait, &target);
        assert!((c - (1.0 - 2.0 * 50.0 / 120.0)).abs() < 1e-9);
    }

    #[test]
    fn test_target_centering_no_overlap_is_zero() {
        let bait = Bait::new("chr1", 100, 200, "b");
        let target = Target::new("chr1", 200, 300, "t");
        assert_eq!(target_centering(&bait, &target), 0.0);
    }

    #[test]
    fn test_target_centering_different_chrom_is_zero() {
        let bait = Bait::new("chr1", 100, 200, "b");
        let target = Target::new("chr2", 100, 200, "t");
        assert_eq!(target_centering(&bait, &target), 0.0);
    }

    #[test]
    fn test_target_centering_snp_at_bait_edge_is_low() {
        // SNP at far left edge of bait [100,220): overlap [100,101),
        // overlap_center=100.5, bait_center=160, bait_len=120 → near 0
        let bait = Bait::new("chr1", 100, 220, "b");
        let target = Target::new("chr1", 100, 101, "snp");
        let c = target_centering(&bait, &target);
        assert!(c < 0.1, "SNP at bait edge should score near 0, got {c}");
    }

    #[test]
    fn test_target_padded() {
        let t = Target::new("chr1", 100, 200, "t");
        let p = t.padded(20);
        assert_eq!(p.start, 80);
        assert_eq!(p.end, 220);
    }

    #[test]
    fn test_target_padded_clamps_to_zero() {
        let t = Target::new("chr1", 5, 50, "t");
        let p = t.padded(100);
        assert_eq!(p.start, 0);
    }

    #[test]
    fn test_parse_bed_error_on_start_greater_than_end() {
        let result = parse_bed_reader(Cursor::new("chr1\t200\t100\tb\n"));
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_interval_list_error_on_zero_start() {
        let input = "@HD\tVN:1.6\nchr1\t0\t10\t+\tbait\n";
        let result = parse_interval_list_reader(Cursor::new(input));
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_interval_list_error_on_start_greater_than_end() {
        // Mirrors the BED parser's start > end guard so inverted coordinates are
        // rejected at parse time instead of underflowing downstream.
        let input = "@HD\tVN:1.6\nchr1\t200\t100\t+\tbad\n";
        let result = parse_interval_list_reader(Cursor::new(input));
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("start (200) > end (100)"),
            "error should report the inverted coordinates"
        );
    }

    #[test]
    fn test_parse_bed_from_path() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t100\t200\tbait").unwrap();
        let baits = parse_bed(f.path()).unwrap();
        assert_eq!(baits.len(), 1);
        assert_eq!(baits[0].name, "bait");
    }

    #[test]
    fn test_parse_interval_list_from_path() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "@HD\tVN:1.6").unwrap();
        writeln!(f, "chr1\t101\t200\t+\tbait").unwrap();
        let baits = parse_interval_list(f.path()).unwrap();
        assert_eq!(baits.len(), 1);
        assert_eq!(baits[0].start, 100);
    }

    #[test]
    fn test_parse_targets_from_path() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t100\t200\ttarget1").unwrap();
        let targets = parse_targets(f.path()).unwrap();
        assert_eq!(targets.len(), 1);
        assert_eq!(targets[0].name, "target1");
    }

    #[test]
    fn test_target_is_empty_true_when_zero_length() {
        assert!(Target::new("chr1", 100, 100, "t").is_empty());
    }

    #[test]
    fn test_target_is_empty_false_when_nonzero_length() {
        assert!(!Target::new("chr1", 100, 101, "t").is_empty());
    }

    #[test]
    fn test_parse_coord_token_start_zero_returns_none() {
        assert!(parse_coord_token("chr1:0-100").is_none());
    }

    #[test]
    fn test_parse_coord_token_start_greater_than_end_returns_none() {
        assert!(parse_coord_token("chr1:200-100").is_none());
    }

    #[test]
    fn test_parse_coord_token_empty_chrom_returns_none() {
        assert!(parse_coord_token(":100-200").is_none());
    }

    #[test]
    fn test_parse_interval_list_error_on_too_few_fields() {
        let input = "@HD\tVN:1.6\nchr1\t100\t200\n";
        assert!(parse_interval_list_reader(Cursor::new(input)).is_err());
    }

    #[test]
    fn test_parse_interval_list_error_on_invalid_start() {
        let input = "@HD\tVN:1.6\nchr1\tabc\t200\t+\tbait\n";
        assert!(parse_interval_list_reader(Cursor::new(input)).is_err());
    }

    #[test]
    fn test_parse_interval_list_error_on_invalid_end() {
        let input = "@HD\tVN:1.6\nchr1\t100\txyz\t+\tbait\n";
        assert!(parse_interval_list_reader(Cursor::new(input)).is_err());
    }

    #[test]
    fn test_parse_bed_error_on_invalid_start() {
        assert!(parse_bed_reader(Cursor::new("chr1\tabc\t100\tb\n")).is_err());
    }

    #[test]
    fn test_parse_bed_error_on_invalid_end() {
        assert!(parse_bed_reader(Cursor::new("chr1\t100\txyz\tb\n")).is_err());
    }

    #[test]
    fn test_parse_interval_list_skips_empty_lines() {
        // An empty line in an Interval List should be silently skipped.
        let input = "@HD\tVN:1.6\n\nchr1\t1\t10\t+\tb\n";
        let baits = parse_interval_list_reader(Cursor::new(input)).unwrap();
        assert_eq!(baits.len(), 1);
        assert_eq!(baits[0].name, "b");
    }

    #[test]
    fn test_parse_interval_list_error_on_four_fields() {
        // A non-header, non-empty line with only 4 fields (missing name) should be an error.
        let input = "chr1\t1\t10\t+\n";
        assert!(parse_interval_list_reader(Cursor::new(input)).is_err());
    }
}
