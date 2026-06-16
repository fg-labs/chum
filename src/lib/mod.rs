//! `chumlib` — hybrid selection bait evaluation library.
#![warn(missing_docs)]

use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use log::*;
use proglog::ProgLogBuilder;
use rayon::prelude::*;

pub mod bait;
pub mod blast;
pub mod intervals;
pub mod mappability;
pub mod metrics;
pub mod repbase;
pub mod rnafold;
pub mod score;
pub mod sequence;

pub use bait::{BaitEvaluator, EvaluatorConfig, write_group_metrics_tsv, write_metrics_tsv};
pub use intervals::{Bait, Target};
pub use metrics::{BaitGroupMetric, BaitMetric};

/// Skip the current test if any of the named tools are absent from `$PATH`.
///
/// Accepts one or more tool names; skips if **any** are missing.
///
/// ```ignore
/// skip_if_missing!("blastn");
/// skip_if_missing!("bgzip", "tabix");
/// ```
#[cfg(test)]
macro_rules! skip_if_missing {
    ($($tool:literal),+ $(,)?) => {
        if $( ::which::which($tool).is_err() )||+ {
            return;
        }
    };
}
#[cfg(test)]
pub(crate) use skip_if_missing;

/// Top-level run arguments, built from the CLI.
#[derive(Debug)]
pub struct ScoreRunArgs {
    /// Input BED, Interval List, or FASTA of baits.
    pub baits: PathBuf,
    /// Output per-bait TSV path (stdout if `None` or `-`).
    pub per_bait: Option<PathBuf>,
    /// Indexed FASTA reference.
    pub reference: Option<PathBuf>,
    /// Target intervals for centering and group metrics.
    pub targets: Option<PathBuf>,
    /// Output per-target group metrics TSV.
    pub per_target: Option<PathBuf>,
    /// Bases to pad targets when matching baits.
    pub target_padding: u32,
    /// BLASTn database name.
    pub blast_db: Option<String>,
    /// Directory containing BLAST databases.
    pub blast_db_path: Option<PathBuf>,
    /// Enable DUST complexity masking.
    pub blast_dust: bool,
    /// Number of CPU threads passed to each `blastn` subprocess via `-num_threads`.
    pub blast_threads: usize,
    /// Mappability bedGraph.
    pub mappability: Option<PathBuf>,
    /// RepBase BED file.
    pub rep_base: Option<PathBuf>,
    /// Enable RNAFold.
    pub oligo_fold: bool,
    /// RNAFold temperature.
    pub oligo_fold_temp: f64,
    /// Bundled ViennaRNA parameter file to use for oligo folding.
    pub oligo_fold_par: rnafold::ParFile,
    /// Parallel worker threads.
    pub threads: usize,
    /// Baits per batch.
    pub batch_size: usize,
}

/// Run `chum` end-to-end given the provided [`ScoreRunArgs`].
pub fn run_score(args: ScoreRunArgs) -> Result<()> {
    // Validate co-required arguments.
    if args.targets.is_some() != args.per_target.is_some() {
        bail!("`--targets` and `--per-target` must both be provided or both omitted");
    }

    if args.batch_size == 0 {
        bail!("`--batch-size` must be greater than 0");
    }

    let optional_files: &[(&str, Option<&PathBuf>)] = &[
        ("--ref", args.reference.as_ref()),
        ("--targets", args.targets.as_ref()),
        ("--mappability", args.mappability.as_ref()),
        ("--rep-base", args.rep_base.as_ref()),
    ];
    for (flag, path) in optional_files
        .iter()
        .filter_map(|(f, p)| p.map(|p| (*f, p)))
    {
        if !path.exists() {
            bail!("{flag} file not found: {}", path.display());
        }
    }

    if !args.baits.exists() {
        bail!("Input file not found: {}", args.baits.display());
    }

    let mut baits = intervals::parse_intervals(&args.baits)
        .with_context(|| format!("Failed to parse input: {}", args.baits.display()))?;

    info!(
        "Read {} bait intervals from {}",
        baits.len(),
        args.baits.display()
    );

    // Load sequences from reference FASTA when provided.
    if let Some(ref fasta_path) = args.reference {
        bait::load_sequences_from_fasta(&mut baits, fasta_path)?;
    } else if baits.iter().any(|b| b.sequence.is_none()) {
        // Only warn when at least one bait lacks a sequence (i.e., input is BED/IL, not FASTA).
        warn!("--ref not provided; sequence-dependent metrics will be skipped");
    }

    // Resolve the RNAFold parameter file from the selected enum variant.
    let (oligo_fold_param, oligo_fold_param_name): (Option<Vec<u8>>, Option<String>) =
        if args.oligo_fold {
            (
                args.oligo_fold_par.bytes().map(|b| b.to_vec()),
                args.oligo_fold_par.file_name().map(|s| s.to_string()),
            )
        } else {
            (None, None)
        };

    let config = EvaluatorConfig {
        reference: args.reference.clone(),
        target_padding: args.target_padding,
        blast_db: args.blast_db.clone(),
        blast_db_path: args.blast_db_path.clone(),
        blast_dust: args.blast_dust,
        blast_threads: args.blast_threads,
        mappability: args.mappability.clone(),
        rep_base: args.rep_base.clone(),
        oligo_fold: args.oligo_fold,
        oligo_fold_temp: args.oligo_fold_temp,
        oligo_fold_param,
        oligo_fold_param_name,
        threads: args.threads,
        batch_size: args.batch_size,
    };

    if args.threads > 1 && args.blast_threads > 1 && args.blast_db.is_some() {
        warn!(
            "BLAST concurrency: {} worker thread(s) × {} blastn thread(s) = up to {} total BLAST threads",
            args.threads,
            args.blast_threads,
            args.threads.saturating_mul(args.blast_threads)
        );
    }

    let evaluator = BaitEvaluator::new(config)?;

    // Fail fast if any located bait sits on a contig the mappability/RepBase indexes
    // don't cover, rather than aborting mid-run on the first affected bait.
    evaluator.validate_bait_contigs(&baits)?;

    let progress = ProgLogBuilder::new()
        .name("main")
        .verb("Evaluated")
        .noun("bait records")
        .unit(100)
        .build();

    let chunks: Vec<&[Bait]> = baits.chunks(args.batch_size).collect();
    let mut all_metrics: Vec<BaitMetric> = Vec::with_capacity(baits.len());

    if args.threads <= 1 {
        for chunk in &chunks {
            let mut chunk_metrics: Vec<BaitMetric> = chunk
                .iter()
                .map(|b| evaluator.evaluate(b))
                .collect::<Result<_>>()?;

            // BLAST is batched per-chunk for efficiency (one subprocess per chunk).
            if let Some(ref runner) = evaluator.blast {
                bait::apply_blast_to_batch(&evaluator, runner, chunk, &mut chunk_metrics)?;
            }

            for _ in chunk.iter() {
                progress.record();
            }
            all_metrics.extend(chunk_metrics);
        }
    } else {
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build()?;
        let results: Vec<Result<Vec<BaitMetric>>> = thread_pool.install(|| {
            chunks
                .par_iter()
                .map(|chunk| {
                    let mut ms: Vec<BaitMetric> = chunk
                        .iter()
                        .map(|b| evaluator.evaluate(b))
                        .collect::<Result<_>>()?;
                    if let Some(ref runner) = evaluator.blast {
                        bait::apply_blast_to_batch(&evaluator, runner, chunk, &mut ms)?;
                    }
                    Ok(ms)
                })
                .collect()
        });
        for result in results {
            let chunk_metrics = result?;
            for _ in &chunk_metrics {
                progress.record();
            }
            all_metrics.extend(chunk_metrics);
        }
    }

    // Parse targets once; used for both centering and group metric output.
    let targets: Vec<Target> = if let Some(ref tpath) = args.targets {
        intervals::parse_targets(tpath)
            .with_context(|| format!("Failed to parse targets: {}", tpath.display()))?
    } else {
        vec![]
    };

    // Apply target-centering metrics.
    // For unlocated FASTA baits, substitute the BLAST top-hit position so that
    // target assignment and centering still work.
    if !targets.is_empty() {
        let located_baits: Vec<Bait> = all_metrics
            .iter()
            .zip(baits.iter())
            .map(|(metric, b)| {
                if b.is_located() {
                    b.clone()
                } else {
                    // Use the pre-computed proxy bait from apply_blast_to_batch when
                    // available, falling back to bait_from_top_hit (parses interval string)
                    // when BLAST was not run or produced no hits.
                    metric
                        .blast_proxy_bait
                        .clone()
                        .or_else(|| bait::bait_from_top_hit(metric))
                        .unwrap_or_else(|| b.clone())
                }
            })
            .collect();
        bait::apply_target_centering(&mut all_metrics, &located_baits, &targets);
    }

    // Write per-bait output.
    let output_writer: Box<dyn std::io::Write> = match &args.per_bait {
        Some(path) if path.to_str().unwrap_or("") != "-" => Box::new(std::fs::File::create(path)?),
        _ => Box::new(std::io::stdout()),
    };
    write_metrics_tsv(output_writer, &all_metrics)?;

    // For unlocated baits, replace the synthetic `unknown:...` interval with the
    // BLAST top-hit interval so group-metric overlap queries work correctly.
    for metric in &mut all_metrics {
        if metric.chrom == "unknown"
            && let Some(proxy) = metric
                .blast_proxy_bait
                .clone()
                .or_else(|| bait::bait_from_top_hit(metric))
        {
            metric.interval = format!("{}:{}-{}", proxy.chrom, proxy.start + 1, proxy.end);
            metric.chrom = proxy.chrom;
            metric.start = proxy.start;
            metric.end = proxy.end;
        }
    }

    if let Some(gpath) = &args.per_target {
        let group_metrics = bait::build_group_metrics(&all_metrics, &targets, args.target_padding);
        info!(
            "Writing {} target group metrics to {}",
            group_metrics.len(),
            gpath.display()
        );
        let group_writer = std::fs::File::create(gpath)?;
        write_group_metrics_tsv(group_writer, &group_metrics)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::Write;
    use std::path::PathBuf;
    use tempfile::NamedTempFile;

    use super::*;

    fn make_args(baits: PathBuf) -> ScoreRunArgs {
        ScoreRunArgs {
            baits,
            per_bait: None,
            reference: None,
            targets: None,
            per_target: None,
            target_padding: 0,
            blast_db: None,
            blast_db_path: None,
            blast_dust: false,
            blast_threads: 1,
            mappability: None,
            rep_base: None,
            oligo_fold: false,
            oligo_fold_temp: 65.0,
            oligo_fold_par: rnafold::ParFile::DnaMathews2004,
            threads: 1,
            batch_size: 50,
        }
    }

    #[test]
    fn test_run_score_nonexistent_baits_errors() {
        let args = make_args(PathBuf::from("/nonexistent/baits.bed"));
        let result = run_score(args);
        assert!(result.is_err());
    }

    #[test]
    fn test_run_score_batch_size_zero_errors() {
        // `--batch-size 0` must be rejected cleanly, not panic in slice::chunks(0).
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t0\t100\tbait1").unwrap();
        let mut args = make_args(f.path().to_path_buf());
        args.batch_size = 0;
        let result = run_score(args);
        assert!(result.is_err());
        assert!(
            result.unwrap_err().to_string().contains("batch-size"),
            "error should mention --batch-size"
        );
    }

    #[test]
    fn test_run_score_targets_without_per_target_errors() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t0\t100\tbait1").unwrap();
        let mut args = make_args(f.path().to_path_buf());
        args.targets = Some(f.path().to_path_buf());
        args.per_target = None;
        let result = run_score(args);
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("--targets"),
            "error should mention --targets, got: {msg}"
        );
    }

    #[test]
    fn test_run_score_per_target_without_targets_errors() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t0\t100\tbait1").unwrap();
        let group_out = NamedTempFile::new().unwrap();
        let mut args = make_args(f.path().to_path_buf());
        args.targets = None;
        args.per_target = Some(group_out.path().to_path_buf());
        let result = run_score(args);
        assert!(result.is_err());
    }

    #[test]
    fn test_run_score_nonexistent_reference_errors() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t0\t100\tbait1").unwrap();
        let mut args = make_args(f.path().to_path_buf());
        args.reference = Some(PathBuf::from("/nonexistent/ref.fa"));
        let result = run_score(args);
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("--ref"),
            "error should mention --ref, got: {msg}"
        );
    }

    #[test]
    fn test_run_score_bed_input_succeeds() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t0\t100\tbait1").unwrap();
        let out = NamedTempFile::new().unwrap();
        let mut args = make_args(f.path().to_path_buf());
        args.per_bait = Some(out.path().to_path_buf());
        let result = run_score(args);
        assert!(
            result.is_ok(),
            "run_score failed: {:?}",
            result.unwrap_err()
        );
        let content = std::fs::read_to_string(out.path()).unwrap();
        assert!(
            content.contains("bait_name"),
            "output missing bait_name header"
        );
        assert!(content.contains("bait1"), "output missing bait name");
    }

    #[test]
    fn test_run_score_outputs_to_stdout_when_no_output_path() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t0\t50\tb1").unwrap();
        let args = make_args(f.path().to_path_buf());
        // output = None → writes to stdout; just verify no error
        let result = run_score(args);
        assert!(
            result.is_ok(),
            "run_score failed: {:?}",
            result.unwrap_err()
        );
    }

    #[test]
    fn test_run_score_parallel_threads_succeeds() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t0\t100\tbait1").unwrap();
        writeln!(f, "chr1\t200\t300\tbait2").unwrap();
        let out = NamedTempFile::new().unwrap();
        let mut args = make_args(f.path().to_path_buf());
        args.threads = 2;
        args.per_bait = Some(out.path().to_path_buf());
        let result = run_score(args);
        assert!(
            result.is_ok(),
            "parallel run_score failed: {:?}",
            result.unwrap_err()
        );
        let content = std::fs::read_to_string(out.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        // header + 2 data rows
        assert_eq!(lines.len(), 3, "expected header + 2 rows");
    }

    #[test]
    fn test_run_score_with_reference_fasta_populates_sequences() {
        // The bundled test FASTA has contig "test-contig"; create matching BED input.
        let fasta = std::path::PathBuf::from("tests/data/blast/hs38DH-chr3:129530791-129531030.fa");
        if !fasta.exists() {
            return;
        }
        let mut bed_f = NamedTempFile::new().unwrap();
        writeln!(bed_f, "test-contig\t0\t30\tbait1").unwrap();
        let out = NamedTempFile::new().unwrap();
        let mut args = make_args(bed_f.path().to_path_buf());
        args.reference = Some(fasta);
        args.per_bait = Some(out.path().to_path_buf());
        let result = run_score(args);
        assert!(
            result.is_ok(),
            "run_score with reference FASTA failed: {:?}",
            result.unwrap_err()
        );
        let content = std::fs::read_to_string(out.path()).unwrap();
        // GC content column should be populated (not empty) since sequence was loaded.
        let lines: Vec<&str> = content.lines().collect();
        let headers: Vec<&str> = lines[0].split('\t').collect();
        let gc_idx = headers.iter().position(|&h| h == "gc_content").unwrap();
        let gc_val = lines[1].split('\t').nth(gc_idx).unwrap();
        assert!(
            !gc_val.is_empty(),
            "gc_content should be non-empty when reference provided"
        );
    }

    #[test]
    fn test_run_score_fasta_without_coords_and_targets() {
        // FASTA bait with no coordinate header → chrom="unknown".
        // Combined with targets, this exercises the "unlocated bait" proxy-bait
        // branch (lines 223-227) and the "fix unknown chrom" loop (lines 244-253).
        let mut fasta_f = NamedTempFile::new().unwrap();
        writeln!(fasta_f, ">bait_no_coords\nACGTACGTACGT").unwrap();
        let mut targets_f = NamedTempFile::new().unwrap();
        writeln!(targets_f, "chr1\t0\t200\ttarget1").unwrap();
        let out = NamedTempFile::new().unwrap();
        let group_out = NamedTempFile::new().unwrap();
        let mut args = make_args(fasta_f.path().to_path_buf());
        args.targets = Some(targets_f.path().to_path_buf());
        args.per_target = Some(group_out.path().to_path_buf());
        args.per_bait = Some(out.path().to_path_buf());
        let result = run_score(args);
        assert!(
            result.is_ok(),
            "run_score with unlocated FASTA bait failed: {:?}",
            result.unwrap_err()
        );
    }

    #[test]
    fn test_run_score_parallel_blast_warning_executes() {
        // With threads > 1, blast_threads > 1, and blast_db set, the warn! fires
        // before BaitEvaluator::new fails (db nonexistent → error is fine).
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t0\t100\tbait1").unwrap();
        let mut args = make_args(f.path().to_path_buf());
        args.threads = 2;
        args.blast_threads = 2;
        args.blast_db = Some("nonexistent_db_for_warning_test".to_string());
        // Just verify the warn! path runs and we get some error (not a panic).
        let _ = run_score(args); // error expected; we only care it doesn't panic
    }

    #[test]
    fn test_run_score_oligo_fold_param_resolved_before_evaluator() {
        // When oligo_fold=true the param resolution lines run even if RNAfold is absent.
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t0\t100\tbait1").unwrap();
        let mut args = make_args(f.path().to_path_buf());
        args.oligo_fold = true;
        args.oligo_fold_par = rnafold::ParFile::DnaMathews2004;
        // If RNAfold is not on PATH, run_score will fail at BaitEvaluator::new —
        // but lines 108-116 (param resolution) still execute first.
        let _ = run_score(args); // may succeed or fail depending on $PATH
    }

    #[test]
    fn test_run_score_with_targets_and_per_target_succeeds() {
        let mut baits_f = NamedTempFile::new().unwrap();
        writeln!(baits_f, "chr1\t0\t100\tbait1").unwrap();
        let mut targets_f = NamedTempFile::new().unwrap();
        writeln!(targets_f, "chr1\t0\t200\ttarget1").unwrap();
        let out = NamedTempFile::new().unwrap();
        let group_out = NamedTempFile::new().unwrap();
        let mut args = make_args(baits_f.path().to_path_buf());
        args.targets = Some(targets_f.path().to_path_buf());
        args.per_target = Some(group_out.path().to_path_buf());
        args.per_bait = Some(out.path().to_path_buf());
        let result = run_score(args);
        assert!(
            result.is_ok(),
            "run_score failed: {:?}",
            result.unwrap_err()
        );
        let group_content = std::fs::read_to_string(group_out.path()).unwrap();
        assert!(
            group_content.contains("target_name"),
            "group output missing header"
        );
        assert!(group_content.contains("target1"));
    }
}
