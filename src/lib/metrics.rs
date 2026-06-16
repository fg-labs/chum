//! Per-bait and per-target aggregated metric types.
use serde::{Serialize, Serializer};
use statrs::statistics::{Data, Max, Min, OrderStatistics};

use crate::blast::BlastHitFormat6;
use crate::intervals::Bait;
use crate::sequence::is_mito_chrom;

/// Serialize an `Option<T>` as an empty string when `None`, or the value when `Some`.
pub fn serialize_option<S, T>(value: &Option<T>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
    T: Serialize,
{
    match value {
        Some(v) => v.serialize(serializer),
        None => serializer.serialize_none(),
    }
}

/// Per-bait quality metrics.
///
/// All optional fields serialize as an empty string when `None`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct BaitMetric {
    /// Name of the bait from the input interval name field.
    pub bait_name: String,

    /// Genomic interval in `chr:start-end` format (1-based closed, matching Picard convention).
    pub interval: String,

    /// Length of the bait sequence in bases.
    pub bait_length: u32,

    /// Reference sequence name (not serialized; use `interval` for display).
    #[serde(skip)]
    pub chrom: String,

    /// 0-based start coordinate (not serialized; use `interval` for display).
    #[serde(skip)]
    pub start: u64,

    /// 0-based end coordinate (not serialized; use `interval` for display).
    #[serde(skip)]
    pub end: u64,

    /// Proxy bait built from the BLAST top-hit coordinates for unlocated FASTA baits.
    /// Populated inside `apply_blast_to_batch`; used for target-centering in lieu of
    /// the original bait when no embedded coordinates were found in the FASTA header.
    #[serde(skip)]
    pub(crate) blast_proxy_bait: Option<crate::intervals::Bait>,

    /// Strand of the bait (`+`, `-`, or `.`), when available from the input or BLAST top hit.
    #[serde(serialize_with = "serialize_option")]
    pub strand: Option<String>,

    /// Fraction of G+C bases in the sequence. `N` and `.` are excluded from the
    /// denominator; lowercase soft-masked bases are case-folded and counted.
    #[serde(serialize_with = "serialize_option")]
    pub gc_content: Option<f64>,

    /// Number of lowercase, N, or `.` bases in the sequence.
    #[serde(serialize_with = "serialize_option")]
    pub masked_bases: Option<u32>,

    /// Count of distinct homopolymer runs of length ≥ 3.
    #[serde(serialize_with = "serialize_option")]
    pub homopolymers_size_3_or_greater: Option<u32>,

    /// Length of the longest homopolymer run.
    #[serde(serialize_with = "serialize_option")]
    pub longest_homopolymer_size: Option<u32>,

    /// Comma-separated RepBase feature names overlapping this bait.
    #[serde(serialize_with = "serialize_option")]
    pub rep_base_features: Option<String>,

    /// Interval of the top BLAST hit in `chr:start-end` format.
    #[serde(serialize_with = "serialize_option")]
    pub blast_top_hit_interval: Option<String>,

    /// `true` if the top BLAST hit interval exactly matches the bait interval.
    #[serde(serialize_with = "serialize_option")]
    pub blast_top_hit_matches: Option<bool>,

    /// `true` if the top BLAST hit interval overlaps the bait interval.
    #[serde(serialize_with = "serialize_option")]
    pub blast_top_hit_overlaps: Option<bool>,

    /// E-value of the top BLAST hit.
    #[serde(serialize_with = "serialize_option")]
    pub blast_top_hit_e_value: Option<f64>,

    /// Percent identity of the top BLAST hit.
    #[serde(serialize_with = "serialize_option")]
    pub blast_top_hit_identity: Option<f64>,

    /// Total number of BLAST hits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_hits: Option<u32>,

    /// Interval of the second-best BLAST hit.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_interval: Option<String>,

    /// Alignment length of the second-best BLAST hit.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_length: Option<u32>,

    /// Percent identity of the second-best BLAST hit.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_percent: Option<f64>,

    /// Gap opens in the second-best BLAST hit.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_gaps: Option<u32>,

    /// E-value of the second-best BLAST hit.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_e_value: Option<f64>,

    /// Interval of the best-scoring BLAST hit to a mitochondrial sequence.
    #[serde(serialize_with = "serialize_option")]
    pub mito_hit_interval: Option<String>,

    /// E-value of the best mitochondrial BLAST hit.
    #[serde(serialize_with = "serialize_option")]
    pub mito_hit_e_value: Option<f64>,

    /// Percent identity of the best mitochondrial BLAST hit.
    #[serde(serialize_with = "serialize_option")]
    pub mito_hit_identity: Option<f64>,

    /// Mean per-base mappability score across the bait.
    #[serde(serialize_with = "serialize_option")]
    pub mean_mappability: Option<f64>,

    /// Minimum per-base mappability score.
    #[serde(serialize_with = "serialize_option")]
    pub min_mappability: Option<f64>,

    /// First-quartile per-base mappability score.
    #[serde(serialize_with = "serialize_option")]
    pub q1_mappability: Option<f64>,

    /// Median per-base mappability score.
    #[serde(serialize_with = "serialize_option")]
    pub median_mappability: Option<f64>,

    /// Third-quartile per-base mappability score.
    #[serde(serialize_with = "serialize_option")]
    pub q3_mappability: Option<f64>,

    /// Maximum per-base mappability score.
    #[serde(serialize_with = "serialize_option")]
    pub max_mappability: Option<f64>,

    /// Number of bases with mappability score equal to 1.0.
    #[serde(serialize_with = "serialize_option")]
    pub unique_mappability: Option<u32>,

    /// Number of bases with mappability score equal to 0.0.
    #[serde(serialize_with = "serialize_option")]
    pub zero_mappability: Option<u32>,

    /// Predicted minimum Gibbs free energy (kcal/mol) from RNAFold.
    #[serde(serialize_with = "serialize_option")]
    pub min_free_energy: Option<f64>,

    /// Predicted secondary structure dot-bracket notation from RNAFold.
    #[serde(serialize_with = "serialize_option")]
    pub folding_structure: Option<String>,

    /// Temperature (°C) used for RNAFold simulation.
    #[serde(serialize_with = "serialize_option")]
    pub oligo_fold_temperature: Option<f64>,

    /// Name of the ViennaRNA parameter file used for RNAFold (e.g. `dna_mathews2004.par`).
    #[serde(serialize_with = "serialize_option")]
    pub oligo_fold_param_file: Option<String>,

    /// Bait nucleotide sequence (uppercase as extracted from reference).
    #[serde(serialize_with = "serialize_option")]
    pub sequence: Option<String>,

    /// Name of the overlapping target interval.
    #[serde(serialize_with = "serialize_option")]
    pub target_name: Option<String>,

    /// Target interval in `chr:start-end` format (1-based closed).
    #[serde(serialize_with = "serialize_option")]
    pub target_interval: Option<String>,

    /// Bait centering over its target: 1.0 = perfectly centered, 0.0 = no overlap.
    #[serde(serialize_with = "serialize_option")]
    pub target_centering: Option<f64>,

    /// Composite quality score from 0 (worst) to 1 (best).
    /// `None` when fewer than two scored components are available.
    #[serde(serialize_with = "serialize_option")]
    pub bait_score: Option<f64>,
}

/// Populate BLAST fields on a [`BaitMetric`] from a slice of format-6 hits.
///
/// The bait's own interval is used to determine whether the top hit "matches"
/// (exact same position) or "overlaps" (any overlap with bait).
pub fn apply_blast_hits(metric: &mut BaitMetric, bait: &Bait, hits: &[BlastHitFormat6]) {
    metric.blast_hits = Some(hits.len() as u32);

    if let Some(top) = hits.first() {
        let top_chrom = &top.sseqid;
        let top_start = (top.sstart.min(top.send) as u64).saturating_sub(1); // to 0-based
        let top_end = top.sstart.max(top.send) as u64;
        let top_interval = format!("{}:{}-{}", top_chrom, top_start + 1, top_end);
        let matches = top_chrom == &bait.chrom && top_start == bait.start && top_end == bait.end;
        let overlaps = top_chrom == &bait.chrom && top_start < bait.end && top_end > bait.start;
        metric.blast_top_hit_interval = Some(top_interval);
        metric.blast_top_hit_matches = Some(matches);
        metric.blast_top_hit_overlaps = Some(overlaps);
        metric.blast_top_hit_e_value = Some(top.evalue);
        metric.blast_top_hit_identity = Some(top.pident);
    }

    // Scan all hits for the best-scoring mitochondrial alignment.
    let best_mito = hits
        .iter()
        .filter(|h| is_mito_chrom(&h.sseqid))
        .min_by(|a, b| {
            a.evalue
                .partial_cmp(&b.evalue)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    if let Some(mito) = best_mito {
        let m_start = (mito.sstart.min(mito.send) as u64).saturating_sub(1);
        let m_end = mito.sstart.max(mito.send) as u64;
        metric.mito_hit_interval = Some(format!("{}:{}-{}", mito.sseqid, m_start + 1, m_end));
        metric.mito_hit_e_value = Some(mito.evalue);
        metric.mito_hit_identity = Some(mito.pident);
    }

    if let Some(second) = hits.get(1) {
        let s_start = (second.sstart.min(second.send) as u64).saturating_sub(1);
        let s_end = second.sstart.max(second.send) as u64;
        metric.blast_second_hit_interval =
            Some(format!("{}:{}-{}", second.sseqid, s_start + 1, s_end));
        metric.blast_second_hit_length = Some(second.length);
        metric.blast_second_hit_percent = Some(second.pident);
        metric.blast_second_hit_gaps = Some(second.gapopen);
        metric.blast_second_hit_e_value = Some(second.evalue);
    }
}

/// Compute mappability summary statistics and populate the relevant fields on `metric`.
///
/// Does nothing when the scores vector is empty.
pub fn apply_mappability(metric: &mut BaitMetric, scores: &[f64]) {
    if scores.is_empty() {
        return;
    }

    let mut data = Data::new(scores.to_vec());

    let mean = scores.iter().sum::<f64>() / scores.len() as f64;
    metric.mean_mappability = Some(mean);
    metric.min_mappability = Some(data.min());
    metric.max_mappability = Some(data.max());
    metric.q1_mappability = Some(data.percentile(25));
    metric.median_mappability = Some(data.percentile(50));
    metric.q3_mappability = Some(data.percentile(75));
    metric.unique_mappability = Some(
        scores
            .iter()
            .filter(|&&s| (s - 1.0).abs() < f64::EPSILON)
            .count() as u32,
    );
    metric.zero_mappability = Some(scores.iter().filter(|&&s| s == 0.0).count() as u32);
}

/// The separator used to join multiple RepBase feature names.
pub const REPBASE_FEATURE_SEPARATOR: &str = ",";

/// Populate the `rep_base_features` field on `metric` from overlapping feature names.
///
/// Sets the field to `None` when there are no overlapping features, or to a
/// comma-joined sorted string otherwise.
pub fn apply_repbase(metric: &mut BaitMetric, features: Vec<String>) {
    if features.is_empty() {
        metric.rep_base_features = None;
    } else {
        metric.rep_base_features = Some(features.join(REPBASE_FEATURE_SEPARATOR));
    }
}

impl BaitMetric {
    /// Construct a metric from required positional fields; all optional fields default to `None`.
    pub fn from_position(bait_name: String, chrom: &str, start: u64, end: u64) -> Self {
        BaitMetric {
            bait_name,
            interval: format!("{}:{}-{}", chrom, start + 1, end),
            bait_length: (end - start) as u32,
            chrom: chrom.to_string(),
            start,
            end,
            blast_proxy_bait: None,
            strand: None,
            gc_content: None,
            masked_bases: None,
            homopolymers_size_3_or_greater: None,
            longest_homopolymer_size: None,
            rep_base_features: None,
            blast_top_hit_interval: None,
            blast_top_hit_matches: None,
            blast_top_hit_overlaps: None,
            blast_top_hit_e_value: None,
            blast_top_hit_identity: None,
            blast_hits: None,
            blast_second_hit_interval: None,
            blast_second_hit_length: None,
            blast_second_hit_percent: None,
            blast_second_hit_gaps: None,
            blast_second_hit_e_value: None,
            mito_hit_interval: None,
            mito_hit_e_value: None,
            mito_hit_identity: None,
            mean_mappability: None,
            min_mappability: None,
            q1_mappability: None,
            median_mappability: None,
            q3_mappability: None,
            max_mappability: None,
            unique_mappability: None,
            zero_mappability: None,
            min_free_energy: None,
            folding_structure: None,
            oligo_fold_temperature: None,
            oligo_fold_param_file: None,
            sequence: None,
            target_name: None,
            target_interval: None,
            target_centering: None,
            bait_score: None,
        }
    }
}

/// Aggregated per-target metrics across all overlapping baits.
///
/// For each numeric `BaitMetric` field the min, mean, and max across overlapping baits are
/// reported. `None` values are excluded from aggregation.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct BaitGroupMetric {
    /// Name of the target interval.
    pub target_name: String,
    /// Target interval in `chr:start-end` format (1-based closed).
    pub target_interval: String,
    /// Length of the target interval in bases.
    pub target_length: u32,
    /// Number of bases of padding applied to both sides of the target.
    pub target_padding: u32,
    /// Number of baits overlapping the (padded) target.
    pub num_baits: u32,

    /// Maximum bait length across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub bait_length_max_of_baits: Option<u32>,
    /// Mean bait length across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub bait_length_mean_of_baits: Option<f64>,
    /// Minimum bait length across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub bait_length_min_of_baits: Option<u32>,

    /// Maximum masked-base count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub masked_bases_max_of_baits: Option<u32>,
    /// Mean masked-base count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub masked_bases_mean_of_baits: Option<f64>,
    /// Minimum masked-base count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub masked_bases_min_of_baits: Option<u32>,

    /// Maximum GC content across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub gc_content_max_of_baits: Option<f64>,
    /// Mean GC content across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub gc_content_mean_of_baits: Option<f64>,
    /// Minimum GC content across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub gc_content_min_of_baits: Option<f64>,

    /// Maximum longest-homopolymer length across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub longest_homopolymer_max_of_baits: Option<u32>,
    /// Mean longest-homopolymer length across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub longest_homopolymer_mean_of_baits: Option<f64>,
    /// Minimum longest-homopolymer length across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub longest_homopolymer_min_of_baits: Option<u32>,

    /// Maximum BLAST hit count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_hits_max_of_baits: Option<u32>,
    /// Mean BLAST hit count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_hits_mean_of_baits: Option<f64>,
    /// Minimum BLAST hit count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_hits_min_of_baits: Option<u32>,

    /// Maximum top-BLAST-hit e-value across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_top_hit_e_value_max_of_baits: Option<f64>,
    /// Mean top-BLAST-hit e-value across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_top_hit_e_value_mean_of_baits: Option<f64>,
    /// Minimum top-BLAST-hit e-value across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_top_hit_e_value_min_of_baits: Option<f64>,

    /// Maximum top-BLAST-hit percent identity across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_top_hit_identity_max_of_baits: Option<f64>,
    /// Mean top-BLAST-hit percent identity across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_top_hit_identity_mean_of_baits: Option<f64>,
    /// Minimum top-BLAST-hit percent identity across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_top_hit_identity_min_of_baits: Option<f64>,

    /// Maximum second-BLAST-hit alignment length across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_length_max_of_baits: Option<u32>,
    /// Mean second-BLAST-hit alignment length across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_length_mean_of_baits: Option<f64>,
    /// Minimum second-BLAST-hit alignment length across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_length_min_of_baits: Option<u32>,

    /// Maximum second-BLAST-hit percent identity across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_percent_max_of_baits: Option<f64>,
    /// Mean second-BLAST-hit percent identity across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_percent_mean_of_baits: Option<f64>,
    /// Minimum second-BLAST-hit percent identity across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_percent_min_of_baits: Option<f64>,

    /// Maximum second-BLAST-hit gap-open count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_gaps_max_of_baits: Option<u32>,
    /// Mean second-BLAST-hit gap-open count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_gaps_mean_of_baits: Option<f64>,
    /// Minimum second-BLAST-hit gap-open count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub blast_second_hit_gaps_min_of_baits: Option<u32>,

    /// Maximum mitochondrial-hit e-value across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub mito_hit_e_value_max_of_baits: Option<f64>,
    /// Mean mitochondrial-hit e-value across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub mito_hit_e_value_mean_of_baits: Option<f64>,
    /// Minimum mitochondrial-hit e-value across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub mito_hit_e_value_min_of_baits: Option<f64>,

    /// Maximum mitochondrial-hit percent identity across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub mito_hit_identity_max_of_baits: Option<f64>,
    /// Mean mitochondrial-hit percent identity across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub mito_hit_identity_mean_of_baits: Option<f64>,
    /// Minimum mitochondrial-hit percent identity across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub mito_hit_identity_min_of_baits: Option<f64>,

    /// Maximum first-quartile mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub q1_mappability_max_of_baits: Option<f64>,
    /// Mean first-quartile mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub q1_mappability_mean_of_baits: Option<f64>,
    /// Minimum first-quartile mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub q1_mappability_min_of_baits: Option<f64>,

    /// Maximum median mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub median_mappability_max_of_baits: Option<f64>,
    /// Mean median mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub median_mappability_mean_of_baits: Option<f64>,
    /// Minimum median mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub median_mappability_min_of_baits: Option<f64>,

    /// Maximum third-quartile mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub q3_mappability_max_of_baits: Option<f64>,
    /// Mean third-quartile mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub q3_mappability_mean_of_baits: Option<f64>,
    /// Minimum third-quartile mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub q3_mappability_min_of_baits: Option<f64>,

    /// Maximum per-bait maximum mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub max_mappability_max_of_baits: Option<f64>,
    /// Mean per-bait maximum mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub max_mappability_mean_of_baits: Option<f64>,
    /// Minimum per-bait maximum mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub max_mappability_min_of_baits: Option<f64>,

    /// Maximum unique-mappability base count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub unique_mappability_max_of_baits: Option<u32>,
    /// Mean unique-mappability base count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub unique_mappability_mean_of_baits: Option<f64>,
    /// Minimum unique-mappability base count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub unique_mappability_min_of_baits: Option<u32>,

    /// Maximum zero-mappability base count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub zero_mappability_max_of_baits: Option<u32>,
    /// Mean zero-mappability base count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub zero_mappability_mean_of_baits: Option<f64>,
    /// Minimum zero-mappability base count across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub zero_mappability_min_of_baits: Option<u32>,

    /// Maximum mean mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub mean_mappability_max_of_baits: Option<f64>,
    /// Mean mean mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub mean_mappability_mean_of_baits: Option<f64>,
    /// Minimum mean mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub mean_mappability_min_of_baits: Option<f64>,

    /// Maximum per-bait minimum mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub min_mappability_max_of_baits: Option<f64>,
    /// Mean per-bait minimum mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub min_mappability_mean_of_baits: Option<f64>,
    /// Minimum per-bait minimum mappability score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub min_mappability_min_of_baits: Option<f64>,

    /// Maximum minimum free energy across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub min_free_energy_max_of_baits: Option<f64>,
    /// Mean minimum free energy across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub min_free_energy_mean_of_baits: Option<f64>,
    /// Minimum minimum free energy across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub min_free_energy_min_of_baits: Option<f64>,

    /// Maximum composite bait score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub bait_score_max_of_baits: Option<f64>,
    /// Mean composite bait score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub bait_score_mean_of_baits: Option<f64>,
    /// Minimum composite bait score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub bait_score_min_of_baits: Option<f64>,

    /// Maximum target-centering score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub target_centering_max_of_baits: Option<f64>,
    /// Mean target-centering score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub target_centering_mean_of_baits: Option<f64>,
    /// Minimum target-centering score across overlapping baits.
    #[serde(serialize_with = "serialize_option")]
    pub target_centering_min_of_baits: Option<f64>,
}

/// Helpers for aggregating slices of optional numeric values.
pub(crate) fn mean_of(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        None
    } else {
        Some(values.iter().sum::<f64>() / values.len() as f64)
    }
}

pub(crate) fn min_of_f64(values: &[f64]) -> Option<f64> {
    values.iter().cloned().reduce(f64::min)
}

pub(crate) fn max_of_f64(values: &[f64]) -> Option<f64> {
    values.iter().cloned().reduce(f64::max)
}

pub(crate) fn min_of_u32(values: &[u32]) -> Option<u32> {
    values.iter().cloned().min()
}

pub(crate) fn max_of_u32(values: &[u32]) -> Option<u32> {
    values.iter().cloned().max()
}

pub(crate) fn mean_of_u32(values: &[u32]) -> Option<f64> {
    if values.is_empty() {
        None
    } else {
        Some(values.iter().map(|v| *v as f64).sum::<f64>() / values.len() as f64)
    }
}

impl BaitGroupMetric {
    /// Aggregate a slice of overlapping [`BaitMetric`] records for a single target.
    pub fn build(
        target_name: String,
        target_interval: String,
        target_length: u32,
        target_padding: u32,
        overlapping: &[&BaitMetric],
    ) -> Self {
        let bait_lengths: Vec<u32> = overlapping.iter().map(|m| m.bait_length).collect();
        let masked: Vec<u32> = overlapping.iter().filter_map(|m| m.masked_bases).collect();
        let gc: Vec<f64> = overlapping.iter().filter_map(|m| m.gc_content).collect();
        let hp: Vec<u32> = overlapping
            .iter()
            .filter_map(|m| m.longest_homopolymer_size)
            .collect();
        let blast_hits: Vec<u32> = overlapping.iter().filter_map(|m| m.blast_hits).collect();
        let b1_eval: Vec<f64> = overlapping
            .iter()
            .filter_map(|m| m.blast_top_hit_e_value)
            .collect();
        let b1_pct: Vec<f64> = overlapping
            .iter()
            .filter_map(|m| m.blast_top_hit_identity)
            .collect();
        let b2_len: Vec<u32> = overlapping
            .iter()
            .filter_map(|m| m.blast_second_hit_length)
            .collect();
        let b2_pct: Vec<f64> = overlapping
            .iter()
            .filter_map(|m| m.blast_second_hit_percent)
            .collect();
        let b2_gaps: Vec<u32> = overlapping
            .iter()
            .filter_map(|m| m.blast_second_hit_gaps)
            .collect();
        let mito_eval: Vec<f64> = overlapping
            .iter()
            .filter_map(|m| m.mito_hit_e_value)
            .collect();
        let mito_pct: Vec<f64> = overlapping
            .iter()
            .filter_map(|m| m.mito_hit_identity)
            .collect();
        let q1_map: Vec<f64> = overlapping
            .iter()
            .filter_map(|m| m.q1_mappability)
            .collect();
        let med_map: Vec<f64> = overlapping
            .iter()
            .filter_map(|m| m.median_mappability)
            .collect();
        let q3_map: Vec<f64> = overlapping
            .iter()
            .filter_map(|m| m.q3_mappability)
            .collect();
        let max_map: Vec<f64> = overlapping
            .iter()
            .filter_map(|m| m.max_mappability)
            .collect();
        let uniq_map: Vec<u32> = overlapping
            .iter()
            .filter_map(|m| m.unique_mappability)
            .collect();
        let zero_map: Vec<u32> = overlapping
            .iter()
            .filter_map(|m| m.zero_mappability)
            .collect();
        let mean_map: Vec<f64> = overlapping
            .iter()
            .filter_map(|m| m.mean_mappability)
            .collect();
        let min_map: Vec<f64> = overlapping
            .iter()
            .filter_map(|m| m.min_mappability)
            .collect();
        let mfe: Vec<f64> = overlapping
            .iter()
            .filter_map(|m| m.min_free_energy)
            .collect();
        let scores: Vec<f64> = overlapping.iter().filter_map(|m| m.bait_score).collect();
        let centering: Vec<f64> = overlapping
            .iter()
            .filter_map(|m| m.target_centering)
            .collect();

        BaitGroupMetric {
            target_name,
            target_interval,
            target_length,
            target_padding,
            num_baits: overlapping.len() as u32,
            bait_length_max_of_baits: max_of_u32(&bait_lengths),
            bait_length_mean_of_baits: mean_of_u32(&bait_lengths),
            bait_length_min_of_baits: min_of_u32(&bait_lengths),
            masked_bases_max_of_baits: max_of_u32(&masked),
            masked_bases_mean_of_baits: mean_of_u32(&masked),
            masked_bases_min_of_baits: min_of_u32(&masked),
            gc_content_max_of_baits: max_of_f64(&gc),
            gc_content_mean_of_baits: mean_of(&gc),
            gc_content_min_of_baits: min_of_f64(&gc),
            longest_homopolymer_max_of_baits: max_of_u32(&hp),
            longest_homopolymer_mean_of_baits: mean_of_u32(&hp),
            longest_homopolymer_min_of_baits: min_of_u32(&hp),
            blast_hits_max_of_baits: max_of_u32(&blast_hits),
            blast_hits_mean_of_baits: mean_of_u32(&blast_hits),
            blast_hits_min_of_baits: min_of_u32(&blast_hits),
            blast_top_hit_e_value_max_of_baits: max_of_f64(&b1_eval),
            blast_top_hit_e_value_mean_of_baits: mean_of(&b1_eval),
            blast_top_hit_e_value_min_of_baits: min_of_f64(&b1_eval),
            blast_top_hit_identity_max_of_baits: max_of_f64(&b1_pct),
            blast_top_hit_identity_mean_of_baits: mean_of(&b1_pct),
            blast_top_hit_identity_min_of_baits: min_of_f64(&b1_pct),
            blast_second_hit_length_max_of_baits: max_of_u32(&b2_len),
            blast_second_hit_length_mean_of_baits: mean_of_u32(&b2_len),
            blast_second_hit_length_min_of_baits: min_of_u32(&b2_len),
            blast_second_hit_percent_max_of_baits: max_of_f64(&b2_pct),
            blast_second_hit_percent_mean_of_baits: mean_of(&b2_pct),
            blast_second_hit_percent_min_of_baits: min_of_f64(&b2_pct),
            blast_second_hit_gaps_max_of_baits: max_of_u32(&b2_gaps),
            blast_second_hit_gaps_mean_of_baits: mean_of_u32(&b2_gaps),
            blast_second_hit_gaps_min_of_baits: min_of_u32(&b2_gaps),
            mito_hit_e_value_max_of_baits: max_of_f64(&mito_eval),
            mito_hit_e_value_mean_of_baits: mean_of(&mito_eval),
            mito_hit_e_value_min_of_baits: min_of_f64(&mito_eval),
            mito_hit_identity_max_of_baits: max_of_f64(&mito_pct),
            mito_hit_identity_mean_of_baits: mean_of(&mito_pct),
            mito_hit_identity_min_of_baits: min_of_f64(&mito_pct),
            q1_mappability_max_of_baits: max_of_f64(&q1_map),
            q1_mappability_mean_of_baits: mean_of(&q1_map),
            q1_mappability_min_of_baits: min_of_f64(&q1_map),
            median_mappability_max_of_baits: max_of_f64(&med_map),
            median_mappability_mean_of_baits: mean_of(&med_map),
            median_mappability_min_of_baits: min_of_f64(&med_map),
            q3_mappability_max_of_baits: max_of_f64(&q3_map),
            q3_mappability_mean_of_baits: mean_of(&q3_map),
            q3_mappability_min_of_baits: min_of_f64(&q3_map),
            max_mappability_max_of_baits: max_of_f64(&max_map),
            max_mappability_mean_of_baits: mean_of(&max_map),
            max_mappability_min_of_baits: min_of_f64(&max_map),
            unique_mappability_max_of_baits: max_of_u32(&uniq_map),
            unique_mappability_mean_of_baits: mean_of_u32(&uniq_map),
            unique_mappability_min_of_baits: min_of_u32(&uniq_map),
            zero_mappability_max_of_baits: max_of_u32(&zero_map),
            zero_mappability_mean_of_baits: mean_of_u32(&zero_map),
            zero_mappability_min_of_baits: min_of_u32(&zero_map),
            mean_mappability_max_of_baits: max_of_f64(&mean_map),
            mean_mappability_mean_of_baits: mean_of(&mean_map),
            mean_mappability_min_of_baits: min_of_f64(&mean_map),
            min_mappability_max_of_baits: max_of_f64(&min_map),
            min_mappability_mean_of_baits: mean_of(&min_map),
            min_mappability_min_of_baits: min_of_f64(&min_map),
            min_free_energy_max_of_baits: max_of_f64(&mfe),
            min_free_energy_mean_of_baits: mean_of(&mfe),
            min_free_energy_min_of_baits: min_of_f64(&mfe),
            bait_score_max_of_baits: max_of_f64(&scores),
            bait_score_mean_of_baits: mean_of(&scores),
            bait_score_min_of_baits: min_of_f64(&scores),
            target_centering_max_of_baits: max_of_f64(&centering),
            target_centering_mean_of_baits: mean_of(&centering),
            target_centering_min_of_baits: min_of_f64(&centering),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bait_metric_from_position() {
        let m = BaitMetric::from_position("bait1".to_string(), "chr1", 100, 200);
        assert_eq!(m.bait_name, "bait1");
        assert_eq!(m.interval, "chr1:101-200");
        assert_eq!(m.bait_length, 100);
        assert!(m.gc_content.is_none());
        assert!(m.sequence.is_none());
    }

    #[test]
    fn test_mean_of_empty_is_none() {
        assert_eq!(mean_of(&[]), None);
        assert_eq!(mean_of_u32(&[]), None);
    }

    #[test]
    fn test_mean_of_values() {
        assert_eq!(mean_of(&[1.0, 3.0]), Some(2.0));
        assert_eq!(mean_of_u32(&[2, 4]), Some(3.0));
    }

    #[test]
    fn test_min_max_of_f64() {
        let v = vec![0.25, 0.75, 0.50];
        assert_eq!(min_of_f64(&v), Some(0.25));
        assert_eq!(max_of_f64(&v), Some(0.75));
    }

    #[test]
    fn test_min_max_of_u32() {
        let v = vec![3u32, 1, 2];
        assert_eq!(min_of_u32(&v), Some(1));
        assert_eq!(max_of_u32(&v), Some(3));
    }

    #[test]
    fn test_bait_group_metric_build_aggregates_correctly() {
        let mut m1 = BaitMetric::from_position("b1".to_string(), "chr1", 10, 30);
        m1.gc_content = Some(0.50);
        m1.longest_homopolymer_size = Some(3);
        let mut m2 = BaitMetric::from_position("b2".to_string(), "chr1", 20, 40);
        m2.gc_content = Some(0.60);
        m2.longest_homopolymer_size = Some(5);

        let refs: Vec<&BaitMetric> = vec![&m1, &m2];
        let g = BaitGroupMetric::build("t1".to_string(), "chr1:11-50".to_string(), 40, 0, &refs);

        assert_eq!(g.num_baits, 2);
        assert_eq!(g.gc_content_min_of_baits, Some(0.50));
        assert_eq!(g.gc_content_max_of_baits, Some(0.60));
        assert_eq!(g.gc_content_mean_of_baits, Some(0.55));
        assert_eq!(g.longest_homopolymer_min_of_baits, Some(3));
        assert_eq!(g.longest_homopolymer_max_of_baits, Some(5));
        assert_eq!(g.longest_homopolymer_mean_of_baits, Some(4.0));
    }

    #[test]
    fn test_bait_group_metric_excludes_none_from_aggregation() {
        let mut m1 = BaitMetric::from_position("b1".to_string(), "chr1", 0, 10);
        m1.gc_content = Some(0.50);
        let mut m2 = BaitMetric::from_position("b2".to_string(), "chr1", 0, 10);
        m2.gc_content = None; // excluded

        let refs: Vec<&BaitMetric> = vec![&m1, &m2];
        let g = BaitGroupMetric::build("t".to_string(), "chr1:1-10".to_string(), 10, 0, &refs);

        assert_eq!(g.num_baits, 2);
        assert_eq!(g.gc_content_min_of_baits, Some(0.50));
        assert_eq!(g.gc_content_max_of_baits, Some(0.50));
        assert_eq!(g.gc_content_mean_of_baits, Some(0.50));
    }

    #[test]
    fn test_apply_blast_hits_no_hits() {
        let bait = Bait::new("chr1", 100, 200, "b");
        let mut metric = BaitMetric::from_position("b".to_string(), "chr1", 100, 200);
        apply_blast_hits(&mut metric, &bait, &[]);
        assert_eq!(metric.blast_hits, Some(0));
        assert!(metric.blast_top_hit_interval.is_none());
        assert!(metric.blast_second_hit_interval.is_none());
    }

    #[test]
    fn test_apply_blast_hits_single_exact_match() {
        // Top hit exactly matches bait chr1:100-200 (0-based) → sstart=101, send=200 (1-based)
        let bait = Bait::new("chr1", 100, 200, "b");
        let mut metric = BaitMetric::from_position("b".to_string(), "chr1", 100, 200);
        let hits = vec![BlastHitFormat6 {
            qseqid: "0".to_string(),
            sseqid: "chr1".to_string(),
            pident: 100.0,
            length: 100,
            mismatch: 0,
            gapopen: 0,
            qstart: 1,
            qend: 100,
            sstart: 101,
            send: 200,
            evalue: 1e-50,
            bitscore: 200.0,
        }];
        apply_blast_hits(&mut metric, &bait, &hits);
        assert_eq!(metric.blast_hits, Some(1));
        assert_eq!(metric.blast_top_hit_matches, Some(true));
        assert_eq!(metric.blast_top_hit_overlaps, Some(true));
        assert!(metric.blast_second_hit_interval.is_none());
    }

    #[test]
    fn test_apply_blast_hits_second_hit_populated() {
        let bait = Bait::new("chr1", 100, 200, "b");
        let mut metric = BaitMetric::from_position("b".to_string(), "chr1", 100, 200);
        let hits = [
            BlastHitFormat6 {
                qseqid: "0".to_string(),
                sseqid: "chr1".to_string(),
                pident: 100.0,
                length: 100,
                mismatch: 0,
                gapopen: 0,
                qstart: 1,
                qend: 100,
                sstart: 101,
                send: 200,
                evalue: 1e-50,
                bitscore: 200.0,
            },
            BlastHitFormat6 {
                qseqid: "0".to_string(),
                sseqid: "chr2".to_string(),
                pident: 95.0,
                length: 98,
                mismatch: 2,
                gapopen: 1,
                qstart: 1,
                qend: 98,
                sstart: 501,
                send: 598,
                evalue: 1e-20,
                bitscore: 150.0,
            },
        ];
        apply_blast_hits(&mut metric, &bait, &hits);
        assert_eq!(metric.blast_hits, Some(2));
        assert_eq!(
            metric.blast_second_hit_interval,
            Some("chr2:501-598".to_string())
        );
        assert_eq!(metric.blast_second_hit_length, Some(98));
        assert!((metric.blast_second_hit_percent.unwrap() - 95.0).abs() < 1e-6);
        assert_eq!(metric.blast_second_hit_gaps, Some(1));
        assert!((metric.blast_second_hit_e_value.unwrap() - 1e-20).abs() < 1e-25);
    }

    #[test]
    fn test_apply_blast_hits_mito_hit_detected() {
        let bait = Bait::new("chr1", 100, 200, "b");
        let mut metric = BaitMetric::from_position("b".to_string(), "chr1", 100, 200);
        let hits = [
            BlastHitFormat6 {
                qseqid: "0".to_string(),
                sseqid: "chr1".to_string(),
                pident: 100.0,
                length: 100,
                mismatch: 0,
                gapopen: 0,
                qstart: 1,
                qend: 100,
                sstart: 101,
                send: 200,
                evalue: 1e-50,
                bitscore: 200.0,
            },
            BlastHitFormat6 {
                qseqid: "0".to_string(),
                sseqid: "chrM".to_string(),
                pident: 90.0,
                length: 95,
                mismatch: 5,
                gapopen: 0,
                qstart: 1,
                qend: 95,
                sstart: 1001,
                send: 1095,
                evalue: 1e-10,
                bitscore: 100.0,
            },
        ];
        apply_blast_hits(&mut metric, &bait, &hits);
        assert_eq!(metric.mito_hit_interval, Some("chrM:1001-1095".to_string()));
        assert!((metric.mito_hit_e_value.unwrap() - 1e-10).abs() < 1e-15);
        assert!((metric.mito_hit_identity.unwrap() - 90.0).abs() < 1e-6);
    }

    #[test]
    fn test_apply_blast_hits_two_mito_hits_best_evalue_wins() {
        // With two mito hits, the min_by comparator is exercised; the hit with
        // the smaller e-value should be selected.
        let bait = Bait::new("chr1", 100, 200, "b");
        let mut metric = BaitMetric::from_position("b".to_string(), "chr1", 100, 200);
        let make_hit = |sseqid: &str, sstart: u32, send: u32, evalue: f64| BlastHitFormat6 {
            qseqid: "0".to_string(),
            sseqid: sseqid.to_string(),
            pident: 90.0,
            length: 100,
            mismatch: 0,
            gapopen: 0,
            qstart: 1,
            qend: 100,
            sstart,
            send,
            evalue,
            bitscore: 100.0,
        };
        let hits = vec![
            make_hit("chr1", 101, 200, 1e-50),    // top hit, non-mito
            make_hit("chrM", 1001, 1100, 1e-5),   // mito hit #1 (worse evalue)
            make_hit("chrMT", 2001, 2100, 1e-20), // mito hit #2 (better evalue)
        ];
        apply_blast_hits(&mut metric, &bait, &hits);
        // The best mito hit should be chrMT with evalue 1e-20.
        assert!(
            metric
                .mito_hit_interval
                .as_deref()
                .unwrap()
                .starts_with("chrMT:"),
            "expected chrMT to win, got {:?}",
            metric.mito_hit_interval
        );
        assert!((metric.mito_hit_e_value.unwrap() - 1e-20).abs() < 1e-25);
    }

    #[test]
    fn test_apply_blast_hits_top_hit_overlaps_not_matches() {
        // Top hit aligns to chr1:105-204 (1-based) = 0-based [104, 204).
        // Bait is [100, 200) → overlaps but coordinates differ → matches=false, overlaps=true.
        let bait = Bait::new("chr1", 100, 200, "b");
        let mut metric = BaitMetric::from_position("b".to_string(), "chr1", 100, 200);
        let hits = [BlastHitFormat6 {
            qseqid: "0".to_string(),
            sseqid: "chr1".to_string(),
            pident: 98.0,
            length: 100,
            mismatch: 2,
            gapopen: 0,
            qstart: 1,
            qend: 100,
            sstart: 105,
            send: 204,
            evalue: 1e-40,
            bitscore: 180.0,
        }];
        apply_blast_hits(&mut metric, &bait, &hits);
        assert_eq!(metric.blast_top_hit_matches, Some(false));
        assert_eq!(metric.blast_top_hit_overlaps, Some(true));
    }

    #[test]
    fn test_apply_blast_hits_top_hit_different_chrom() {
        // Top hit is on a different chromosome → neither matches nor overlaps.
        let bait = Bait::new("chr1", 100, 200, "b");
        let mut metric = BaitMetric::from_position("b".to_string(), "chr1", 100, 200);
        let hits = [BlastHitFormat6 {
            qseqid: "0".to_string(),
            sseqid: "chr2".to_string(),
            pident: 100.0,
            length: 100,
            mismatch: 0,
            gapopen: 0,
            qstart: 1,
            qend: 100,
            sstart: 101,
            send: 200,
            evalue: 1e-50,
            bitscore: 200.0,
        }];
        apply_blast_hits(&mut metric, &bait, &hits);
        assert_eq!(metric.blast_top_hit_matches, Some(false));
        assert_eq!(metric.blast_top_hit_overlaps, Some(false));
    }

    #[test]
    fn test_apply_mappability_empty_does_nothing() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 10);
        apply_mappability(&mut m, &[]);
        assert!(m.mean_mappability.is_none());
    }

    #[test]
    fn test_apply_mappability_known_values() {
        // scores 0.00, 0.25, 0.75, 1.00 → 4 bases
        // mean=0.50, min=0.00, max=1.00, unique=1, zero=1
        let scores = vec![0.00, 0.25, 0.75, 1.00];
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 4);
        apply_mappability(&mut m, &scores);
        assert!((m.mean_mappability.unwrap() - 0.50).abs() < 1e-9);
        assert!((m.min_mappability.unwrap() - 0.00).abs() < 1e-9);
        assert!((m.max_mappability.unwrap() - 1.00).abs() < 1e-9);
        assert_eq!(m.unique_mappability, Some(1));
        assert_eq!(m.zero_mappability, Some(1));
    }

    #[test]
    fn test_apply_mappability_all_unique() {
        let scores = vec![1.0; 10];
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 10);
        apply_mappability(&mut m, &scores);
        assert_eq!(m.unique_mappability, Some(10));
        assert_eq!(m.zero_mappability, Some(0));
        assert!((m.mean_mappability.unwrap() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_apply_repbase_empty_sets_none() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 10);
        apply_repbase(&mut m, vec![]);
        assert!(m.rep_base_features.is_none());
    }

    #[test]
    fn test_apply_repbase_single_feature() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 10);
        apply_repbase(&mut m, vec!["Alu".to_string()]);
        assert_eq!(m.rep_base_features, Some("Alu".to_string()));
    }

    #[test]
    fn test_apply_repbase_multiple_features_joined() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 10);
        apply_repbase(&mut m, vec!["(TAACCC)n".to_string(), "L1MC5a".to_string()]);
        assert_eq!(m.rep_base_features, Some("(TAACCC)n,L1MC5a".to_string()));
    }

    #[test]
    fn test_serialize_option_as_empty_for_none() {
        #[derive(Serialize)]
        struct Row {
            #[serde(serialize_with = "serialize_option")]
            val: Option<f64>,
        }
        let row = Row { val: None };
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(vec![]);
        wtr.serialize(&row).unwrap();
        let out = String::from_utf8(wtr.into_inner().unwrap()).unwrap();
        // The csv crate quotes empty fields as `""` (QuoteStyle::Necessary).
        assert_eq!(out.trim(), "\"\"");
    }

    #[test]
    fn test_serialize_option_as_value_for_some() {
        #[derive(Serialize)]
        struct Row {
            #[serde(serialize_with = "serialize_option")]
            val: Option<f64>,
        }
        let row = Row { val: Some(0.5) };
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(vec![]);
        wtr.serialize(&row).unwrap();
        let out = String::from_utf8(wtr.into_inner().unwrap()).unwrap();
        assert_eq!(out.trim(), "0.5");
    }
}
