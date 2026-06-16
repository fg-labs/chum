//! Bait evaluator — pre-loaded evaluation state and batch processing helpers.

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::PathBuf;
use std::sync::Mutex;

use anyhow::{Context, Result, bail};
use log::*;
use rust_lapper::{Interval, Lapper};

use crate::intervals::{Bait, Target};
use crate::metrics::{BaitGroupMetric, BaitMetric};
use crate::{blast, intervals, mappability, metrics, repbase, rnafold, score, sequence};

/// Configuration for constructing a [`BaitEvaluator`].
#[derive(Debug, Clone)]
pub struct EvaluatorConfig {
    /// Indexed FASTA reference for sequence-based metrics.
    pub reference: Option<PathBuf>,
    /// Number of bases to pad targets when collecting group metrics.
    pub target_padding: u32,
    /// BLASTn database name.
    pub blast_db: Option<String>,
    /// Directory containing BLAST databases.
    pub blast_db_path: Option<PathBuf>,
    /// Enable DUST complexity masking in BLAST.
    pub blast_dust: bool,
    /// Number of CPU threads passed to each `blastn` subprocess via `-num_threads`.
    pub blast_threads: usize,
    /// Block-compressed bedGraph of mappability scores.
    pub mappability: Option<PathBuf>,
    /// Block-compressed BED of RepBase features.
    pub rep_base: Option<PathBuf>,
    /// Enable RNAFold secondary structure prediction.
    pub oligo_fold: bool,
    /// Temperature in °C for RNAFold.
    pub oligo_fold_temp: f64,
    /// ViennaRNA parameter file content; `None` = RNAfold built-in RNA defaults.
    /// Defaults to the vendored `dna_mathews2004.par` bytes.
    pub oligo_fold_param: Option<Vec<u8>>,
    /// Display name of the parameter file written to the output TSV (e.g. `"dna_mathews2004.par"`).
    /// `None` when RNAFold is disabled.
    pub oligo_fold_param_name: Option<String>,
    /// Number of parallel worker threads.
    pub threads: usize,
    /// Number of baits per batch.
    pub batch_size: usize,
}

/// Pre-loaded state for evaluating baits.
///
/// Expensive resources (BLAST configuration, mappability reader, RepBase reader,
/// optional RNAFold process) are loaded once at construction time and then
/// amortized across all bait evaluations. `BaitEvaluator` is designed to be
/// constructed once and reused from an external bait-designer loop:
///
/// ```no_run
/// use chumlib::{BaitEvaluator, EvaluatorConfig, Bait};
///
/// let config = EvaluatorConfig {
///     reference: None,
///     target_padding: 0,
///     blast_db: None,
///     blast_db_path: None,
///     blast_dust: false,
///     blast_threads: 1,
///     mappability: None,
///     rep_base: None,
///     oligo_fold: false,
///     oligo_fold_temp: 65.0,
///     oligo_fold_param: None,
///     oligo_fold_param_name: None,
///     threads: 1,
///     batch_size: 50,
/// };
/// let evaluator = BaitEvaluator::new(config).unwrap();
/// let bait = Bait::new("chr1", 100, 220, "my-bait");
/// let metric = evaluator.evaluate(&bait).unwrap();
/// ```
pub struct BaitEvaluator {
    pub(crate) config: EvaluatorConfig,
    pub(crate) blast: Option<blast::BlastRunner>,
    pub(crate) mappability: Option<mappability::MappabilityReader>,
    pub(crate) repbase: Option<repbase::RepBaseReader>,
    pub(crate) rnafold: Option<Mutex<rnafold::RnaFoldProcess>>,
}

impl std::fmt::Debug for BaitEvaluator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BaitEvaluator")
            .field("config", &self.config)
            .field("blast", &self.blast)
            .field(
                "mappability_path",
                &self.mappability.as_ref().map(|_| "<MappabilityReader>"),
            )
            .field(
                "repbase_path",
                &self.repbase.as_ref().map(|_| "<RepBaseReader>"),
            )
            .field(
                "rnafold",
                &self.rnafold.as_ref().map(|_| "<RnaFoldProcess>"),
            )
            .finish()
    }
}

impl BaitEvaluator {
    /// Construct and pre-load all resources from `config`.
    pub fn new(config: EvaluatorConfig) -> Result<Self> {
        // Validate that blastn is available when a database is requested.
        if let Some(ref db) = config.blast_db {
            if which::which("blastn").is_err() {
                bail!("`blastn` is not on $PATH but --blast-db was provided");
            }
            if !blast::database_exists(db, config.blast_db_path.as_deref()) {
                bail!(
                    "BLAST database `{}` does not exist{}",
                    db,
                    config
                        .blast_db_path
                        .as_ref()
                        .map(|p| format!(", or in the path at `{}`", p.display()))
                        .unwrap_or_default()
                );
            }
        }

        // Validate that RNAfold is available when requested.
        if config.oligo_fold && which::which("RNAfold").is_err() {
            bail!("`RNAfold` is not on $PATH but --oligo-fold was requested");
        }

        let blast = config.blast_db.as_ref().map(|db| {
            blast::BlastRunner::new(
                db,
                config.blast_db_path.clone(),
                config.blast_dust,
                config.blast_threads,
            )
        });

        let mappability = match &config.mappability {
            Some(path) => Some(
                mappability::MappabilityReader::new(path)
                    .with_context(|| format!("Cannot load mappability: {}", path.display()))?,
            ),
            None => None,
        };

        let repbase = match &config.rep_base {
            Some(path) => Some(
                repbase::RepBaseReader::new(path)
                    .with_context(|| format!("Cannot load RepBase: {}", path.display()))?,
            ),
            None => None,
        };

        let rnafold = if config.oligo_fold {
            let param = config.oligo_fold_param.as_deref();
            Some(Mutex::new(rnafold::RnaFoldProcess::spawn(
                config.oligo_fold_temp,
                param,
            )?))
        } else {
            None
        };

        Ok(BaitEvaluator {
            config,
            blast,
            mappability,
            repbase,
            rnafold,
        })
    }

    /// Evaluate a single bait and return its metrics (without BLAST).
    ///
    /// BLAST is applied separately via [`apply_blast_to_batch`] so that a single
    /// subprocess call handles a whole chunk of baits efficiently.
    pub fn evaluate(&self, bait: &Bait) -> Result<BaitMetric> {
        let mut metric =
            BaitMetric::from_position(bait.name.clone(), &bait.chrom, bait.start, bait.end);
        metric.strand = bait.strand.map(|c| c.to_string());

        // Sequence metrics require a sequence stored on the bait.
        if let Some(seq) = &bait.sequence {
            metric.gc_content = Some(sequence::gc_content(seq));
            metric.masked_bases = Some(sequence::masked_bases(seq));
            metric.homopolymers_size_3_or_greater =
                Some(sequence::homopolymers_size_3_or_greater(seq));
            metric.longest_homopolymer_size = Some(sequence::longest_homopolymer_size(seq));
            metric.sequence = Some(seq.to_ascii_uppercase());
        } else if self.config.reference.is_some() {
            warn!(
                "Bait '{}' has no sequence loaded; sequence metrics skipped",
                bait.name
            );
        }

        // RepBase — only when the bait has real genomic coordinates.
        if let Some(rb) = &self.repbase
            && bait.is_located()
        {
            let features = rb.overlapping_features(bait)?;
            metrics::apply_repbase(&mut metric, features);
        }

        // Mappability — only when the bait has real genomic coordinates.
        if let Some(mi) = &self.mappability
            && bait.is_located()
        {
            let scores = mi.scores_for_bait(bait)?;
            metrics::apply_mappability(&mut metric, &scores);
        }

        // Setup RNAFold.
        if let Some(folder) = &self.rnafold
            && let Some(seq) = &bait.sequence
        {
            let result = folder.lock().unwrap().fold(seq)?;
            metric.min_free_energy = Some(result.min_free_energy);
            metric.folding_structure = Some(result.structure);
            metric.oligo_fold_temperature = Some(self.config.oligo_fold_temp);
            metric.oligo_fold_param_file = self.config.oligo_fold_param_name.clone();
        }

        // Preliminary composite score (without BLAST; recomputed after blast batch).
        metric.bait_score = score::bait_score(&metric);

        Ok(metric)
    }

    /// Validate that every located bait's contig is present in the configured
    /// mappability and RepBase tabix indexes, failing fast with an actionable message.
    ///
    /// This catches a missing contig or a naming mismatch (e.g. `chr1` vs `1`) up front,
    /// before any baits are evaluated, instead of aborting mid-run on the first affected
    /// bait with a low-level tabix error. Unlocated (FASTA) baits are skipped here; their
    /// coordinates are resolved later from the BLAST top hit.
    pub(crate) fn validate_bait_contigs(&self, baits: &[Bait]) -> Result<()> {
        if self.mappability.is_none() && self.repbase.is_none() {
            return Ok(());
        }
        let mut contigs: Vec<&str> = baits
            .iter()
            .filter(|b| b.is_located())
            .map(|b| b.chrom.as_str())
            .collect();
        contigs.sort_unstable();
        contigs.dedup();
        if contigs.is_empty() {
            return Ok(());
        }

        let missing_in = |names: &HashSet<String>| -> Vec<String> {
            contigs
                .iter()
                .copied()
                .filter(|c| !names.contains(*c))
                .map(|c| c.to_string())
                .collect()
        };

        if let Some(mi) = &self.mappability {
            let missing = missing_in(&mi.reference_names()?);
            if !missing.is_empty() {
                bail!(
                    "bait contig(s) not present in the mappability index{}: {}. \
                     Check that contig names match the index (e.g. `chr1` vs `1`) and that \
                     the index covers every bait contig.",
                    self.config
                        .mappability
                        .as_ref()
                        .map(|p| format!(" ({})", p.display()))
                        .unwrap_or_default(),
                    missing.join(", "),
                );
            }
        }

        if let Some(rb) = &self.repbase {
            let missing = missing_in(&rb.reference_names()?);
            if !missing.is_empty() {
                bail!(
                    "bait contig(s) not present in the RepBase index{}: {}. \
                     Check that contig names match the index (e.g. `chr1` vs `1`) and that \
                     the index covers every bait contig.",
                    self.config
                        .rep_base
                        .as_ref()
                        .map(|p| format!(" ({})", p.display()))
                        .unwrap_or_default(),
                    missing.join(", "),
                );
            }
        }

        Ok(())
    }
}

/// Populate BLAST fields for a batch of baits in one subprocess call and recompute scores.
///
/// For FASTA-input baits without embedded coordinates (`!bait.is_located()`), the BLAST
/// top-hit interval is parsed back into a temporary [`Bait`] and used to query mappability
/// and RepBase — so those metrics still populate even when coordinates weren't in the header.
pub(crate) fn apply_blast_to_batch(
    evaluator: &BaitEvaluator,
    runner: &blast::BlastRunner,
    baits: &[Bait],
    metric_batch: &mut [BaitMetric],
) -> Result<()> {
    let hit_batches = runner.align_batch(baits)?;
    for (metric, (bait, hits)) in metric_batch
        .iter_mut()
        .zip(baits.iter().zip(hit_batches.iter()))
    {
        metrics::apply_blast_hits(metric, bait, hits);

        // Derive strand from the BLAST top hit when the input didn't supply one.
        if metric.strand.is_none()
            && let Some(top) = hits.first()
        {
            metric.strand = Some(if top.sstart <= top.send {
                "+".to_string()
            } else {
                "-".to_string()
            });
        }

        // For unlocated baits, retroactively apply positional metrics using the top hit.
        if !bait.is_located()
            && let Some(proxy) = bait_from_top_hit(metric)
        {
            // Consistent with the located-bait path and up-front validation: a missing
            // contig is raised, not silently swallowed.
            if let Some(ref rb) = evaluator.repbase {
                let features = rb.overlapping_features(&proxy).with_context(|| {
                    format!(
                        "RepBase lookup failed for BLAST top-hit contig '{}'",
                        proxy.chrom
                    )
                })?;
                metrics::apply_repbase(metric, features);
            }
            if let Some(ref mi) = evaluator.mappability {
                let scores = mi.scores_for_bait(&proxy).with_context(|| {
                    format!(
                        "mappability lookup failed for BLAST top-hit contig '{}'",
                        proxy.chrom
                    )
                })?;
                metrics::apply_mappability(metric, &scores);
            }
            // Store the proxy so mod.rs can use it for target-centering without
            // re-parsing the blast_top_hit_interval string.
            metric.blast_proxy_bait = Some(proxy);
        }

        metric.bait_score = score::bait_score(metric);
    }
    Ok(())
}

/// Build a temporary [`Bait`] from the BLAST top-hit interval stored on a metric.
///
/// The interval is stored as `chrom:start-end` in 1-based closed (Picard) convention.
/// Returns `None` if the field is absent or the string cannot be parsed.
pub(crate) fn bait_from_top_hit(metric: &BaitMetric) -> Option<Bait> {
    let interval = metric.blast_top_hit_interval.as_deref()?;
    let (chrom, range) = interval.split_once(':')?;
    let (start1, end1) = range.split_once('-')?;
    let start: u64 = start1.parse::<u64>().ok()?.checked_sub(1)?; // 1-based → 0-based
    let end: u64 = end1.parse().ok()?;
    if start > end {
        return None;
    }
    Some(Bait::new(chrom, start, end, &metric.bait_name))
}

/// Load bait sequences from an indexed FASTA file using `bio::io::fasta::IndexedReader`.
///
/// Coordinates are 0-based half-open (BED convention), matching what `bio` expects for `fetch`.
pub(crate) fn load_sequences_from_fasta(
    baits: &mut [Bait],
    fasta_path: &std::path::Path,
) -> Result<()> {
    use bio::io::fasta::IndexedReader;
    use std::path::PathBuf;

    let fai_path = {
        let mut p = fasta_path.as_os_str().to_owned();
        p.push(".fai");
        PathBuf::from(p)
    };
    if !fai_path.exists() {
        bail!(
            "FASTA index not found: {}; run `samtools faidx {}` to create it",
            fai_path.display(),
            fasta_path.display()
        );
    }

    let mut reader = IndexedReader::from_file(&fasta_path)
        .with_context(|| format!("Cannot open indexed FASTA: {}", fasta_path.display()))?;

    for bait in baits.iter_mut() {
        if bait.sequence.is_some() {
            continue; // sequence already loaded from FASTA input — do not re-fetch
        }
        match reader.fetch(&bait.chrom, bait.start, bait.end) {
            Ok(()) => {
                let mut seq_bytes = Vec::new();
                if let Err(e) = reader.read(&mut seq_bytes) {
                    warn!("Cannot read FASTA sequence for bait '{}': {e}", bait.name);
                    continue;
                }
                let seq = String::from_utf8_lossy(&seq_bytes).into_owned();
                // Probes on the minus strand hybridize to the forward strand;
                // the probe sequence itself is the reverse complement.
                bait.sequence = Some(if bait.strand == Some('-') {
                    sequence::reverse_complement(&seq).with_context(|| {
                        format!(
                            "Cannot reverse-complement sequence for bait '{}'",
                            bait.name
                        )
                    })?
                } else {
                    seq
                });
            }
            Err(e) => {
                warn!("Cannot fetch FASTA region for bait '{}': {e}", bait.name);
            }
        }
    }
    Ok(())
}

/// Populate target-centering fields on metrics using a lapper tree.
pub(crate) fn apply_target_centering(
    metrics: &mut [BaitMetric],
    baits: &[Bait],
    targets: &[Target],
) {
    // Build per-chromosome lappers from targets.
    let mut by_chrom: HashMap<String, Vec<Interval<u64, usize>>> = HashMap::new();
    for (idx, target) in targets.iter().enumerate() {
        by_chrom
            .entry(target.chrom.clone())
            .or_default()
            .push(Interval {
                start: target.start,
                stop: target.end,
                val: idx,
            });
    }
    let lappers: HashMap<String, Lapper<u64, usize>> = by_chrom
        .into_iter()
        .map(|(c, ivs)| (c, Lapper::new(ivs)))
        .collect();

    for (metric, bait) in metrics.iter_mut().zip(baits.iter()) {
        let Some(lapper) = lappers.get(&bait.chrom) else {
            continue;
        };
        let hits: Vec<&Interval<u64, usize>> = lapper.find(bait.start, bait.end).collect();
        if hits.is_empty() {
            continue;
        }
        let mut names = Vec::with_capacity(hits.len());
        let mut target_intervals = Vec::with_capacity(hits.len());
        let mut best_centering: f64 = 0.0;
        for hit in &hits {
            let target = &targets[hit.val];
            names.push(target.name.clone());
            target_intervals.push(format!(
                "{}:{}-{}",
                target.chrom,
                target.start + 1,
                target.end
            ));
            let c = intervals::target_centering(bait, target);
            if c > best_centering {
                best_centering = c;
            }
        }
        metric.target_name = Some(names.join(","));
        metric.target_interval = Some(target_intervals.join(","));
        metric.target_centering = Some(best_centering);
    }
}

/// Build group metrics from bait metrics and target intervals.
pub(crate) fn build_group_metrics(
    metrics: &[BaitMetric],
    targets: &[Target],
    padding: u32,
) -> Vec<BaitGroupMetric> {
    // Build a per-chrom lapper from bait metric intervals.
    let mut by_chrom: HashMap<String, Vec<Interval<u64, usize>>> = HashMap::new();
    for (i, m) in metrics.iter().enumerate() {
        by_chrom.entry(m.chrom.clone()).or_default().push(Interval {
            start: m.start,
            stop: m.end,
            val: i,
        });
    }
    let lappers: HashMap<String, Lapper<u64, usize>> = by_chrom
        .into_iter()
        .map(|(c, ivs)| (c, Lapper::new(ivs)))
        .collect();

    let mut seen_baits: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let group_metrics: Vec<BaitGroupMetric> = targets
        .iter()
        .map(|target| {
            let padded = target.padded(padding);
            let lapper = lappers.get(&padded.chrom);
            let overlapping_ivs: Vec<usize> = lapper
                .map(|l| l.find(padded.start, padded.end).map(|iv| iv.val).collect())
                .unwrap_or_default();
            let overlapping_refs: Vec<&BaitMetric> =
                overlapping_ivs.iter().map(|&i| &metrics[i]).collect();

            if overlapping_refs.is_empty() {
                warn!(
                    "Target {} (padding = {}bp) has no overlapping baits.",
                    target.name, padding
                );
            }

            for &i in &overlapping_ivs {
                seen_baits.insert(i);
            }

            let target_interval = format!("{}:{}-{}", target.chrom, target.start + 1, target.end);
            BaitGroupMetric::build(
                target.name.clone(),
                target_interval,
                target.len() as u32,
                padding,
                &overlapping_refs,
            )
        })
        .collect();

    if seen_baits.is_empty() {
        warn!("No baits overlapped with any targets.");
    }
    for (i, m) in metrics.iter().enumerate() {
        if !seen_baits.contains(&i) {
            warn!(
                "Bait with name '{}' was not assigned to any target.",
                m.bait_name
            );
        }
    }

    group_metrics
}

/// Write a slice of serializable records to a tab-separated writer.
fn write_tsv<W: Write, T: serde::Serialize>(writer: W, records: &[T]) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(writer);
    for r in records {
        wtr.serialize(r)?;
    }
    wtr.flush()?;
    Ok(())
}

/// Write a slice of [`BaitMetric`] records to a tab-separated writer.
///
/// The CSV header is written automatically from struct field names on the first
/// `serialize()` call.
pub fn write_metrics_tsv<W: Write>(writer: W, metrics: &[BaitMetric]) -> Result<()> {
    write_tsv(writer, metrics)
}

/// Write a slice of [`BaitGroupMetric`] records to a tab-separated writer.
pub fn write_group_metrics_tsv<W: Write>(writer: W, metrics: &[BaitGroupMetric]) -> Result<()> {
    write_tsv(writer, metrics)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::skip_if_missing;

    fn minimal_config() -> EvaluatorConfig {
        EvaluatorConfig {
            reference: None,
            target_padding: 0,
            blast_db: None,
            blast_db_path: None,
            blast_dust: false,
            blast_threads: 1,
            mappability: None,
            rep_base: None,
            oligo_fold: false,
            oligo_fold_temp: 65.0,
            oligo_fold_param: None,
            oligo_fold_param_name: None,
            threads: 1,
            batch_size: 50,
        }
    }

    #[test]
    fn test_bait_evaluator_debug_contains_struct_name() {
        let evaluator = BaitEvaluator::new(minimal_config()).unwrap();
        let s = format!("{evaluator:?}");
        assert!(s.contains("BaitEvaluator"));
        assert!(s.contains("config"));
    }

    #[test]
    fn test_bait_from_top_hit_valid_interval() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr2", 0, 100);
        m.blast_top_hit_interval = Some("chr2:101-200".to_string());
        let bait = bait_from_top_hit(&m).unwrap();
        assert_eq!(bait.chrom, "chr2");
        assert_eq!(bait.start, 100);
        assert_eq!(bait.end, 200);
    }

    #[test]
    fn test_bait_from_top_hit_no_interval_returns_none() {
        let m = BaitMetric::from_position("b".to_string(), "chr1", 0, 100);
        assert!(bait_from_top_hit(&m).is_none());
    }

    #[test]
    fn test_bait_from_top_hit_unparseable_string_returns_none() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 100);
        m.blast_top_hit_interval = Some("not-a-valid-interval".to_string());
        assert!(bait_from_top_hit(&m).is_none());
    }

    #[test]
    fn test_bait_from_top_hit_start_greater_than_end_returns_none() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 100);
        m.blast_top_hit_interval = Some("chr1:200-100".to_string());
        assert!(bait_from_top_hit(&m).is_none());
    }

    #[test]
    fn test_bait_from_top_hit_zero_start_returns_none() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 100);
        m.blast_top_hit_interval = Some("chr1:0-100".to_string());
        assert!(bait_from_top_hit(&m).is_none());
    }

    #[test]
    fn test_apply_target_centering_assigns_target_name_and_centering() {
        let mut metrics = vec![BaitMetric::from_position(
            "b1".to_string(),
            "chr1",
            150,
            270,
        )];
        let baits = vec![Bait::new("chr1", 150, 270, "b1")];
        let targets = vec![Target::new("chr1", 100, 300, "target1")];
        apply_target_centering(&mut metrics, &baits, &targets);
        assert_eq!(metrics[0].target_name.as_deref(), Some("target1"));
        assert!(metrics[0].target_centering.is_some());
        assert!((metrics[0].target_centering.unwrap() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_apply_target_centering_different_chrom_leaves_fields_empty() {
        let mut metrics = vec![BaitMetric::from_position("b".to_string(), "chr2", 0, 100)];
        let baits = vec![Bait::new("chr2", 0, 100, "b")];
        let targets = vec![Target::new("chr1", 0, 100, "t")];
        apply_target_centering(&mut metrics, &baits, &targets);
        assert!(metrics[0].target_name.is_none());
        assert!(metrics[0].target_centering.is_none());
    }

    #[test]
    fn test_apply_target_centering_multiple_targets_uses_best_centering() {
        let mut metrics = vec![BaitMetric::from_position("b".to_string(), "chr1", 100, 200)];
        let baits = vec![Bait::new("chr1", 100, 200, "b")];
        let targets = vec![
            Target::new("chr1", 80, 220, "big"),
            Target::new("chr1", 110, 130, "small"),
        ];
        apply_target_centering(&mut metrics, &baits, &targets);
        assert!(metrics[0].target_name.is_some());
        assert_eq!(metrics[0].target_centering, Some(1.0));
    }

    #[test]
    fn test_build_group_metrics_no_bait_overlap_produces_zero_bait_count() {
        let metrics = vec![BaitMetric::from_position("b".to_string(), "chr1", 500, 600)];
        let targets = vec![Target::new("chr1", 0, 100, "lonely_target")];
        let result = build_group_metrics(&metrics, &targets, 0);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].num_baits, 0);
    }

    #[test]
    fn test_build_group_metrics_correct_overlap_counts_bait() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 50, 150);
        m.gc_content = Some(0.50);
        let targets = vec![Target::new("chr1", 100, 200, "t1")];
        let result = build_group_metrics(&[m], &targets, 0);
        assert_eq!(result[0].num_baits, 1);
        assert_eq!(result[0].gc_content_max_of_baits, Some(0.50));
    }

    #[test]
    fn test_build_group_metrics_padding_extends_overlap_window() {
        let m = BaitMetric::from_position("b".to_string(), "chr1", 0, 50);
        let targets = vec![Target::new("chr1", 100, 200, "t1")];
        let result_no_pad = build_group_metrics(std::slice::from_ref(&m), &targets, 0);
        let result_padded = build_group_metrics(&[m], &targets, 60);
        assert_eq!(result_no_pad[0].num_baits, 0);
        assert_eq!(result_padded[0].num_baits, 1);
    }

    #[test]
    fn test_evaluate_with_sequence_populates_sequence_metrics() {
        let evaluator = BaitEvaluator::new(minimal_config()).unwrap();
        let bait = Bait::with_sequence("chr1", 0, 12, "b1", "GCGCACGTTTTT");
        let metric = evaluator.evaluate(&bait).unwrap();
        assert!(metric.gc_content.is_some());
        assert!(metric.masked_bases.is_some());
        assert!(metric.homopolymers_size_3_or_greater.is_some());
        assert!(metric.longest_homopolymer_size.is_some());
        assert_eq!(metric.sequence.as_deref(), Some("GCGCACGTTTTT"));
    }

    #[test]
    fn test_evaluate_without_sequence_leaves_sequence_metrics_empty() {
        let evaluator = BaitEvaluator::new(minimal_config()).unwrap();
        let bait = Bait::new("chr1", 0, 100, "b1");
        let metric = evaluator.evaluate(&bait).unwrap();
        assert!(metric.gc_content.is_none());
        assert!(metric.masked_bases.is_none());
        assert!(metric.homopolymers_size_3_or_greater.is_none());
        assert!(metric.longest_homopolymer_size.is_none());
    }

    #[test]
    fn test_evaluate_propagates_strand_from_bait() {
        let evaluator = BaitEvaluator::new(minimal_config()).unwrap();
        let mut bait = Bait::new("chr1", 0, 100, "b1");
        bait.strand = Some('+');
        let metric = evaluator.evaluate(&bait).unwrap();
        assert_eq!(metric.strand.as_deref(), Some("+"));
    }

    #[test]
    fn test_evaluate_sequence_is_uppercased() {
        let evaluator = BaitEvaluator::new(minimal_config()).unwrap();
        let bait = Bait::with_sequence("chr1", 0, 6, "b", "acgtnn");
        let metric = evaluator.evaluate(&bait).unwrap();
        assert_eq!(metric.sequence.as_deref(), Some("ACGTNN"));
    }

    #[test]
    fn test_write_metrics_tsv_produces_header_and_row() {
        let metrics = vec![BaitMetric::from_position("b1".to_string(), "chr1", 0, 100)];
        let mut output = Vec::new();
        write_metrics_tsv(&mut output, &metrics).unwrap();
        let s = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert!(lines.len() >= 2, "expected header + at least one row");
        assert!(lines[0].contains("bait_name"), "header missing bait_name");
        assert!(lines[1].contains("b1"), "data row missing bait name");
    }

    #[test]
    fn test_write_group_metrics_tsv_produces_header_and_row() {
        let m = BaitMetric::from_position("b1".to_string(), "chr1", 50, 150);
        let targets = vec![Target::new("chr1", 0, 200, "t1")];
        let groups = build_group_metrics(&[m], &targets, 0);
        let mut output = Vec::new();
        write_group_metrics_tsv(&mut output, &groups).unwrap();
        let s = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert!(lines.len() >= 2, "expected header + at least one row");
        assert!(
            lines[0].contains("target_name"),
            "header missing target_name"
        );
        assert!(lines[1].contains("t1"), "data row missing target name");
    }

    #[test]
    fn test_evaluate_reference_set_but_no_sequence_still_succeeds() {
        // When a reference is configured but the bait has no sequence, evaluate()
        // emits a warning and returns a metric with None sequence fields (no panic).
        let mut config = minimal_config();
        config.reference = Some(std::path::PathBuf::from("/fake/reference.fa"));
        let evaluator = BaitEvaluator::new(config).unwrap();
        let bait = Bait::new("chr1", 0, 100, "b_no_seq");
        let metric = evaluator.evaluate(&bait).unwrap();
        assert!(metric.gc_content.is_none());
        assert!(metric.sequence.is_none());
    }

    #[test]
    fn test_bait_evaluator_new_errors_when_rnafold_not_on_path() {
        if which::which("RNAfold").is_ok() {
            return; // skip when RNAfold is installed
        }
        let mut config = minimal_config();
        config.oligo_fold = true;
        let result = BaitEvaluator::new(config);
        assert!(result.is_err());
        assert!(
            result.unwrap_err().to_string().contains("RNAfold"),
            "error should mention RNAfold"
        );
    }

    #[test]
    fn test_bait_evaluator_new_errors_when_blastn_not_on_path() {
        if which::which("blastn").is_ok() {
            return; // skip when blastn is installed
        }
        let mut config = minimal_config();
        config.blast_db = Some("some_db".to_string());
        let result = BaitEvaluator::new(config);
        assert!(result.is_err());
        assert!(
            result.unwrap_err().to_string().contains("blastn"),
            "error should mention blastn"
        );
    }

    #[test]
    fn test_bait_evaluator_new_errors_when_blast_db_missing_but_blastn_on_path() {
        skip_if_missing!("blastn");
        let mut config = minimal_config();
        config.blast_db = Some("definitely_nonexistent_blast_db_xyz123".to_string());
        let result = BaitEvaluator::new(config);
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("does not exist"),
            "error should mention 'does not exist', got: {msg}"
        );
    }

    #[test]
    fn test_bait_from_top_hit_single_base_interval() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 1);
        m.blast_top_hit_interval = Some("chrX:50-50".to_string());
        let bait = bait_from_top_hit(&m).unwrap();
        assert_eq!(bait.chrom, "chrX");
        assert_eq!(bait.start, 49); // 1-based → 0-based
        assert_eq!(bait.end, 50);
    }

    #[test]
    fn test_apply_target_centering_assigns_interval_string() {
        let mut metrics = vec![BaitMetric::from_position("b".to_string(), "chr1", 10, 20)];
        let baits = vec![Bait::new("chr1", 10, 20, "b")];
        let targets = vec![Target::new("chr1", 0, 30, "t1")];
        apply_target_centering(&mut metrics, &baits, &targets);
        let interval = metrics[0].target_interval.as_deref().unwrap();
        // Interval List format: chr:start+1-end (1-based closed)
        assert!(
            interval.starts_with("chr1:"),
            "interval should start with chr1:"
        );
    }

    #[test]
    fn test_build_group_metrics_multiple_targets() {
        let m1 = BaitMetric::from_position("b1".to_string(), "chr1", 10, 50);
        let m2 = BaitMetric::from_position("b2".to_string(), "chr1", 110, 150);
        let targets = vec![
            Target::new("chr1", 0, 100, "t1"),
            Target::new("chr1", 100, 200, "t2"),
        ];
        let groups = build_group_metrics(&[m1, m2], &targets, 0);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].num_baits, 1);
        assert_eq!(groups[1].num_baits, 1);
    }

    #[test]
    fn test_load_sequences_from_fasta_error_on_missing_fai() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, ">chr1\nACGT").unwrap();
        let mut baits = vec![Bait::new("chr1", 0, 4, "b")];
        let result = load_sequences_from_fasta(&mut baits, f.path());
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("FASTA index not found")
        );
    }

    /// Path to the pre-built indexed reference FASTA bundled with test data.
    fn blast_test_fasta() -> std::path::PathBuf {
        std::path::PathBuf::from("tests/data/blast/hs38DH-chr3:129530791-129531030.fa")
    }

    #[test]
    fn test_load_sequences_from_fasta_populates_sequence() {
        let fasta = blast_test_fasta();
        if !fasta.exists() {
            return;
        }
        let mut baits = vec![Bait::new("test-contig", 0, 30, "b")];
        load_sequences_from_fasta(&mut baits, &fasta).unwrap();
        let seq = baits[0]
            .sequence
            .as_deref()
            .expect("sequence should be set");
        assert_eq!(seq.len(), 30, "expected 30 bases");
        // First 30 bases of the test-contig FASTA.
        assert_eq!(&seq[..4], "CTAG");
    }

    #[test]
    fn test_load_sequences_from_fasta_skips_bait_with_existing_sequence() {
        let fasta = blast_test_fasta();
        if !fasta.exists() {
            return;
        }
        // Bait already has a sequence; load_sequences_from_fasta must not overwrite it.
        let mut baits = vec![Bait::with_sequence("test-contig", 0, 4, "b", "AAAA")];
        load_sequences_from_fasta(&mut baits, &fasta).unwrap();
        assert_eq!(baits[0].sequence.as_deref(), Some("AAAA"));
    }

    #[test]
    fn test_load_sequences_from_fasta_minus_strand_reverse_complements() {
        let fasta = blast_test_fasta();
        if !fasta.exists() {
            return;
        }
        let mut fwd = vec![Bait::new("test-contig", 0, 10, "fwd")];
        let mut rev = vec![{
            let mut b = Bait::new("test-contig", 0, 10, "rev");
            b.strand = Some('-');
            b
        }];
        load_sequences_from_fasta(&mut fwd, &fasta).unwrap();
        load_sequences_from_fasta(&mut rev, &fasta).unwrap();
        let forward_seq = fwd[0].sequence.as_deref().expect("forward seq");
        let reverse_seq = rev[0].sequence.as_deref().expect("reverse seq");
        assert_eq!(
            reverse_seq,
            crate::sequence::reverse_complement(forward_seq).unwrap()
        );
    }

    #[test]
    fn test_load_sequences_from_fasta_unknown_chrom_leaves_sequence_none() {
        let fasta = blast_test_fasta();
        if !fasta.exists() {
            return;
        }
        // Nonexistent chromosome — fetch fails, warn is emitted, sequence stays None.
        let mut baits = vec![Bait::new("nonexistent_chrom_xyz", 0, 10, "b")];
        load_sequences_from_fasta(&mut baits, &fasta).unwrap();
        assert!(baits[0].sequence.is_none());
    }

    fn mappability_fixture() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mappability.bedgraph.gz")
    }

    fn repbase_fixture() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/repbase.bed.gz")
    }

    #[test]
    fn test_bait_evaluator_new_with_mappability() {
        let mut config = minimal_config();
        config.mappability = Some(mappability_fixture());
        let result = BaitEvaluator::new(config);
        assert!(
            result.is_ok(),
            "BaitEvaluator::new with mappability failed: {:?}",
            result.unwrap_err()
        );
    }

    #[test]
    fn test_validate_bait_contigs_missing_contig_errors() {
        // The mappability fixture covers only chr1; a bait on chr2 must fail validation.
        let mut config = minimal_config();
        config.mappability = Some(mappability_fixture());
        let evaluator = BaitEvaluator::new(config).unwrap();
        let baits = vec![Bait::new("chr2", 0, 100, "b")];
        let result = evaluator.validate_bait_contigs(&baits);
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("chr2"),
            "error should name the missing contig, got: {msg}"
        );
        assert!(
            msg.contains("mappability"),
            "error should mention the mappability index, got: {msg}"
        );
    }

    #[test]
    fn test_validate_bait_contigs_present_contig_ok() {
        let mut config = minimal_config();
        config.mappability = Some(mappability_fixture());
        let evaluator = BaitEvaluator::new(config).unwrap();
        let baits = vec![Bait::new("chr1", 0, 100, "b")];
        assert!(evaluator.validate_bait_contigs(&baits).is_ok());
    }

    #[test]
    fn test_validate_bait_contigs_skips_unlocated_baits() {
        // Unlocated (FASTA) baits have chrom "unknown" and are resolved later from BLAST;
        // they must not trip up-front validation.
        let mut config = minimal_config();
        config.mappability = Some(mappability_fixture());
        let evaluator = BaitEvaluator::new(config).unwrap();
        let baits = vec![Bait::new("unknown", 0, 100, "b")];
        assert!(evaluator.validate_bait_contigs(&baits).is_ok());
    }

    #[test]
    fn test_validate_bait_contigs_no_index_is_ok() {
        let evaluator = BaitEvaluator::new(minimal_config()).unwrap();
        let baits = vec![Bait::new("whatever", 0, 100, "b")];
        assert!(evaluator.validate_bait_contigs(&baits).is_ok());
    }

    #[test]
    fn test_bait_evaluator_new_with_repbase() {
        let mut config = minimal_config();
        config.rep_base = Some(repbase_fixture());
        let result = BaitEvaluator::new(config);
        assert!(
            result.is_ok(),
            "BaitEvaluator::new with repbase failed: {:?}",
            result.unwrap_err()
        );
    }

    #[test]
    fn test_evaluate_with_repbase_populates_features() {
        // Fixture: chr1 0-50 Alu, chr1 0-100 L1, chr1 50-100 Alu, chr1 500-600 MIR3.
        // Bait [50,150) overlaps L1 [0,100) and Alu [50,100) → "Alu,L1".
        let mut config = minimal_config();
        config.rep_base = Some(repbase_fixture());
        let evaluator = BaitEvaluator::new(config).unwrap();
        let bait = Bait::new("chr1", 50, 150, "b");
        let metric = evaluator.evaluate(&bait).unwrap();
        assert_eq!(
            metric.rep_base_features.as_deref(),
            Some("Alu,L1"),
            "RepBase features should be populated"
        );
    }

    #[test]
    fn test_evaluate_with_mappability_populates_scores() {
        // Fixture: chr1 0-200 0.75. Bait [50,150) is fully covered → mean 0.75.
        let mut config = minimal_config();
        config.mappability = Some(mappability_fixture());
        let evaluator = BaitEvaluator::new(config).unwrap();
        let bait = Bait::new("chr1", 50, 150, "b");
        let metric = evaluator.evaluate(&bait).unwrap();
        assert!(
            metric.mean_mappability.is_some(),
            "mean_mappability should be populated"
        );
        assert!(
            (metric.mean_mappability.unwrap() - 0.75).abs() < 1e-9,
            "mean_mappability should be 0.75"
        );
    }

    #[test]
    fn test_bait_evaluator_debug_with_mappability_and_repbase() {
        let mut config = minimal_config();
        config.mappability = Some(mappability_fixture());
        config.rep_base = Some(repbase_fixture());
        let evaluator = BaitEvaluator::new(config).unwrap();
        let debug_str = format!("{evaluator:?}");
        assert!(
            debug_str.contains("MappabilityReader"),
            "debug output should show MappabilityReader"
        );
        assert!(
            debug_str.contains("RepBaseReader"),
            "debug output should show RepBaseReader"
        );
    }

    #[test]
    fn test_bait_evaluator_new_with_rnafold_succeeds_if_available() {
        skip_if_missing!("RNAfold");
        let mut config = minimal_config();
        config.oligo_fold = true;
        config.oligo_fold_param = Some(crate::rnafold::DNA_MATHEWS2004_PAR.to_vec());
        let result = BaitEvaluator::new(config);
        assert!(
            result.is_ok(),
            "BaitEvaluator::new with RNAfold failed: {:?}",
            result.unwrap_err()
        );
    }

    #[test]
    fn test_evaluate_with_rnafold_populates_mfe_if_available() {
        skip_if_missing!("RNAfold");
        let mut config = minimal_config();
        config.oligo_fold = true;
        config.oligo_fold_temp = 65.0;
        config.oligo_fold_param = Some(crate::rnafold::DNA_MATHEWS2004_PAR.to_vec());
        config.oligo_fold_param_name = Some("dna_mathews2004.par".to_string());
        let evaluator = BaitEvaluator::new(config).unwrap();
        let bait = Bait::with_sequence("chr1", 0, 16, "b", "GCGCGCTTTTGCGCGC");
        let metric = evaluator.evaluate(&bait).unwrap();
        assert!(metric.min_free_energy.is_some(), "MFE should be populated");
        assert!(
            metric.folding_structure.is_some(),
            "structure should be populated"
        );
        assert_eq!(metric.oligo_fold_temperature, Some(65.0));
        assert_eq!(
            metric.oligo_fold_param_file.as_deref(),
            Some("dna_mathews2004.par")
        );
    }

    #[test]
    fn test_bait_evaluator_debug_with_rnafold_some_if_available() {
        skip_if_missing!("RNAfold");
        let mut config = minimal_config();
        config.oligo_fold = true;
        config.oligo_fold_param = Some(crate::rnafold::DNA_MATHEWS2004_PAR.to_vec());
        let evaluator = BaitEvaluator::new(config).unwrap();
        let debug_str = format!("{evaluator:?}");
        assert!(
            debug_str.contains("RnaFoldProcess"),
            "debug output should mention RnaFoldProcess"
        );
    }

    /// Build a [`BlastRunner`] pointing at the pre-built BLAST database bundled with test data.
    fn blast_runner_fixture() -> blast::BlastRunner {
        let db_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/blast");
        blast::BlastRunner::new("hs38DH-chr3:129530791-129531030", Some(db_path), false, 1)
    }

    /// First 60 bp of the `test-contig` sequence in the bundled BLAST test FASTA.
    const TEST_SEQ_60: &str = "CTAGCTACCCTCTCCCTGTCTAGGGGGGAGTGCACCCTCCTTAGGCAGTGGGGTCTGTGC";

    #[test]
    fn test_apply_blast_to_batch_empty_batch_returns_ok() {
        // Empty slices must succeed without spawning blastn.
        let evaluator = BaitEvaluator::new(minimal_config()).unwrap();
        let runner = blast::BlastRunner::new("dummy", None, false, 1);
        let result = apply_blast_to_batch(&evaluator, &runner, &[], &mut []);
        assert!(result.is_ok());
    }

    #[test]
    fn test_apply_blast_to_batch_populates_blast_fields() {
        skip_if_missing!("blastn");
        let runner = blast_runner_fixture();
        let evaluator = BaitEvaluator::new(minimal_config()).unwrap();
        let seq = TEST_SEQ_60.to_string();
        let bait = Bait::with_sequence("test-contig", 0, seq.len() as u64, "b0", seq);
        let mut metrics = vec![evaluator.evaluate(&bait).unwrap()];
        apply_blast_to_batch(&evaluator, &runner, &[bait], &mut metrics).unwrap();
        assert!(
            metrics[0].blast_hits.is_some(),
            "blast_hits should be populated after apply_blast_to_batch"
        );
        assert!(
            metrics[0].blast_top_hit_interval.is_some(),
            "blast_top_hit_interval should be populated after apply_blast_to_batch"
        );
    }

    #[test]
    fn test_apply_blast_to_batch_derives_forward_strand_from_hit() {
        skip_if_missing!("blastn");
        let runner = blast_runner_fixture();
        let evaluator = BaitEvaluator::new(minimal_config()).unwrap();
        let seq = TEST_SEQ_60.to_string();
        let bait = Bait::with_sequence("test-contig", 0, seq.len() as u64, "b0", seq);
        let mut metric = evaluator.evaluate(&bait).unwrap();
        metric.strand = None; // ensure unset before BLAST
        let mut metrics = vec![metric];
        apply_blast_to_batch(&evaluator, &runner, &[bait], &mut metrics).unwrap();
        assert_eq!(
            metrics[0].strand.as_deref(),
            Some("+"),
            "forward BLAST hit should set strand to '+'"
        );
    }

    #[test]
    fn test_apply_blast_to_batch_derives_reverse_strand_from_hit() {
        skip_if_missing!("blastn");
        let runner = blast_runner_fixture();
        let evaluator = BaitEvaluator::new(minimal_config()).unwrap();
        let rc_seq = crate::sequence::reverse_complement(TEST_SEQ_60).unwrap();
        let len = rc_seq.len() as u64;
        // RC of the forward sequence should align on the minus strand (sstart > send).
        let bait = Bait::with_sequence("test-contig", 0, len, "b0", rc_seq);
        let mut metric = evaluator.evaluate(&bait).unwrap();
        metric.strand = None;
        let mut metrics = vec![metric];
        apply_blast_to_batch(&evaluator, &runner, &[bait], &mut metrics).unwrap();
        assert_eq!(
            metrics[0].strand.as_deref(),
            Some("-"),
            "reverse-complement BLAST hit should set strand to '-'"
        );
    }

    #[test]
    fn test_apply_blast_to_batch_preserves_existing_strand() {
        skip_if_missing!("blastn");
        let runner = blast_runner_fixture();
        let evaluator = BaitEvaluator::new(minimal_config()).unwrap();
        let seq = TEST_SEQ_60.to_string();
        let mut bait = Bait::with_sequence("test-contig", 0, seq.len() as u64, "b0", seq);
        bait.strand = Some('+');
        let mut metrics = vec![evaluator.evaluate(&bait).unwrap()];
        apply_blast_to_batch(&evaluator, &runner, &[bait], &mut metrics).unwrap();
        assert_eq!(
            metrics[0].strand.as_deref(),
            Some("+"),
            "pre-set strand should not be overwritten by BLAST top hit"
        );
    }

    #[test]
    fn test_apply_blast_to_batch_unlocated_bait_sets_proxy_bait() {
        skip_if_missing!("blastn");
        let runner = blast_runner_fixture();
        let evaluator = BaitEvaluator::new(minimal_config()).unwrap();
        let seq = TEST_SEQ_60.to_string();
        let len = seq.len() as u64;
        // chrom = "unknown" → is_located() == false
        let bait = Bait::with_sequence("unknown", 0, len, "b0", seq);
        let mut metrics = vec![evaluator.evaluate(&bait).unwrap()];
        apply_blast_to_batch(&evaluator, &runner, &[bait], &mut metrics).unwrap();
        assert!(
            metrics[0].blast_proxy_bait.is_some(),
            "unlocated bait with BLAST hits should have a proxy bait stored"
        );
        let proxy = metrics[0].blast_proxy_bait.as_ref().unwrap();
        assert_eq!(
            proxy.chrom, "test-contig",
            "proxy bait chrom should match BLAST top-hit subject"
        );
    }
}
