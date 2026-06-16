//! Composite bait quality scoring.
//!
//! The `bait_score` function is a three-stage product:
//!
//! ```text
//! bait_score = quality_score × blast_multiplier × mito_multiplier
//! ```
//!
//! **Quality score** is a weighted average of synthesis/structural metrics:
//!
//! | Component           | Weight | Formula                                             |
//! |---------------------|-------:|-----------------------------------------------------|
//! | GC content          |      2 | 1.0 if 40–60 %, 0.5 if 30–40/60–70 %, else 0.0      |
//! | Longest homopolymer |      1 | `max(0, 1 − (len − 2) / 8)`                         |
//! | Secondary structure |      1 | 1.0 if MFE ≥ −9; quadratic decay to 0.0 at −15      |
//! | Mean mappability    |      2 | `mean_mappability` directly                         |
//!
//! **BLAST multiplier** is applied as a gate that can drive the final score to zero:
//!
//! | Sub-component              | Formula                                             |
//! |----------------------------|-----------------------------------------------------|
//! | Hit-count factor           | `clamp(1 − log10(hits) / 2, 0, 1)` (0 at ≥100 hits) |
//! | Second-hit identity factor | `1 − clamp((pct − 75) / 25, 0, 1)` (0 at 100 %).    |
//! | Combined                   | `hit_factor × identity_factor`                      |
//!
//! When BLAST is not run the multiplier is 1.0 (no penalty). A score of `None`
//! is returned when fewer than two quality components are available and BLAST has
//! not been run, since a single data point provides little comparative value.
//!
//! **Mito multiplier** penalizes a high-identity off-target hit to a mitochondrial
//! sequence (NUMTs are a common false-positive source):
//!
//! | Condition                                       | Factor                                          |
//! |-------------------------------------------------|-------------------------------------------------|
//! | Bait or its target is on a mito chromosome      | 1.0 (no penalty)                                |
//! | Best mito hit is also the second-best BLAST hit | 1.0 (already penalized by `blast_multiplier`)   |
//! | Otherwise, from `mito_hit_identity`             | `1 − clamp((pct − 75) / 25, 0, 1)` (0 at 100 %) |
//!
//! The second condition prevents double-counting the same alignment's identity in
//! both the BLAST and mito multipliers.
//!
//! ## Relationship to IDT probe scoring
//!
//! IDT's OligoAnalyzer complexity score (scale 1–500) reflects synthesis difficulty
//! for RNA-capture probes. Higher IDT scores indicate probes that pass more of IDT's
//! internal quality filters (GC content, hairpin avoidance, repeat masking, etc.).
//! Probes with an IDT score of 500 are considered "perfect" by IDT's metric, while
//! probes scored at 1 carry the most synthesis risk.
//!
//! The `bait_score` below is a *complementary* specificity/mappability score based on
//! the genomic features computed here. Use both together: an IDT score ensures
//! synthesis quality while `bait_score` captures genomic uniqueness and structural
//! simplicity.

use crate::metrics::BaitMetric;
use crate::sequence::is_mito_chrom;

/// Compute a composite quality score for `metric` in the range [0, 1].
///
/// Returns `None` when fewer than two scored components are available, since a
/// single-component score provides little comparative value.
///
/// # Examples
/// ```
/// use chumlib::{metrics::BaitMetric, score::bait_score};
///
/// let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 120);
/// m.gc_content = Some(0.50);
/// m.longest_homopolymer_size = Some(3);
/// let s = bait_score(&m);
/// assert!(s.is_some());
/// assert!(*s.as_ref().unwrap() > 0.0);
/// assert!(*s.as_ref().unwrap() <= 1.0);
/// ```
pub fn bait_score(metric: &BaitMetric) -> Option<f64> {
    let mut weight_total = 0.0_f64;
    let mut weighted = 0.0_f64;
    let mut components = 0u32;

    // ── GC content (weight 2) ────────────────────────────────────────────────
    if let Some(gc) = metric.gc_content {
        weighted += 2.0 * gc_score(gc);
        weight_total += 2.0;
        components += 1;
    }

    // ── Longest homopolymer (weight 1) ───────────────────────────────────────
    if let Some(hp) = metric.longest_homopolymer_size {
        weighted += 1.0 * homopolymer_score(hp);
        weight_total += 1.0;
        components += 1;
    }

    // ── Secondary structure / MFE (weight 1) ────────────────────────────────
    if let Some(mfe) = metric.min_free_energy {
        weighted += 1.0 * mfe_score(mfe);
        weight_total += 1.0;
        components += 1;
    }

    // ── Mappability (weight 2) ───────────────────────────────────────────────
    if let Some(mm) = metric.mean_mappability {
        weighted += 2.0 * mm.clamp(0.0, 1.0);
        weight_total += 2.0;
        components += 1;
    }

    // Need at least two quality components, or one quality component + BLAST,
    // to produce a score worth comparing.
    let has_blast = metric.blast_hits.is_some();
    if components == 0 || (components < 2 && !has_blast) || weight_total == 0.0 {
        return None;
    }

    let quality_score = weighted / weight_total;

    // ── BLAST multiplier (gates the quality score) ───────────────────────────
    // A non-specific bait drives the final score toward zero regardless of how
    // good its synthesis metrics are.
    let blast_multiplier = match metric.blast_hits {
        None => 1.0,
        Some(hits) => blast_specificity_score(hits, metric.blast_second_hit_percent),
    };

    // ── Mito multiplier ──────────────────────────────────────────────────────
    // No penalty when the bait itself is known to be on a mito chromosome
    // (located BED/Interval List input whose chrom is e.g. chrM).
    let bait_on_mito = is_mito_chrom(&metric.chrom);

    // No penalty when the bait's assigned target(s) are on mito — covers FASTA
    // baits that lack embedded coordinates but overlap a known mito target.
    let target_on_mito = metric
        .target_interval
        .as_deref()
        .map(|ti| {
            ti.split(',').any(|iv| {
                iv.split_once(':')
                    .map(|(c, _)| is_mito_chrom(c))
                    .unwrap_or(false)
            })
        })
        .unwrap_or(false);

    // Avoid double-counting the identity penalty: when the best mito hit IS the
    // second-best BLAST hit, its identity was already penalized by `blast_multiplier`
    // (via `blast_second_hit_percent`). Applying the mito penalty again would square
    // the same factor, so skip it in that case.
    let mito_is_second_hit = metric.mito_hit_interval.is_some()
        && metric.mito_hit_interval == metric.blast_second_hit_interval;

    let mito_multiplier = if bait_on_mito || target_on_mito || mito_is_second_hit {
        1.0
    } else {
        match metric.mito_hit_identity {
            Some(pct) => 1.0 - ((pct - 75.0) / 25.0).clamp(0.0, 1.0),
            None => 1.0,
        }
    };

    Some(quality_score * blast_multiplier * mito_multiplier)
}

/// GC-content component score.
///
/// Returns 1.0 for optimal (40–60 %), 0.5 for marginal (30–40 or 60–70 %), and 0.0
/// for extreme values outside 30–70 %.
pub fn gc_score(gc: f64) -> f64 {
    if (0.40..=0.60).contains(&gc) {
        1.0
    } else if (0.30..0.40).contains(&gc) || (0.60..0.70).contains(&gc) {
        0.5
    } else {
        0.0
    }
}

/// Homopolymer component score based on the longest run length.
///
/// A run of ≤ 2 scores 1.0; the score decreases linearly to 0.0 at run length 10.
pub fn homopolymer_score(longest: u32) -> f64 {
    (1.0 - (longest as f64 - 2.0).max(0.0) / 8.0).max(0.0)
}

/// Secondary-structure component score from the minimum free energy (kcal/mol).
///
/// Structures with MFE ≥ −9 kcal/mol score 1.0 (no meaningful secondary structure).
/// Below −9 the score falls quadratically to 0.0 at −15 kcal/mol, producing a steep
/// penalty for stable hairpins. Values below −15 clamp to 0.0.
///
/// | MFE (kcal/mol) | Score |
/// |---------------:|------:|
/// | ≥ −9           | 1.000 |
/// | −12            | 0.250 |
/// | ≤ −15          | 0.000 |
pub fn mfe_score(mfe: f64) -> f64 {
    if mfe >= -9.0 {
        1.0
    } else if mfe <= -15.0 {
        0.0
    } else {
        // Quadratic fall: t=0 at -9 (score 1.0), t=1 at -15 (score 0.0).
        let t = (mfe + 9.0) / (-6.0);
        (1.0 - t) * (1.0 - t)
    }
}

/// BLAST-specificity multiplier combining hit-count and second-hit identity.
///
/// **Hit-count factor** (log-scale): 1.0 at 1 hit, 0.5 at 10 hits, 0.0 at ≥ 100 hits.
///
/// **Second-hit identity factor**: penalizes a high-identity off-target hit.
/// No penalty below 75 % identity; linearly falls to 0.0 at 100 % identity.
/// When there is no second hit the factor is 1.0 (no penalty).
///
/// The two factors are multiplied, so a bait with 74,000 hits and a 92 %
/// second hit scores 0.0 × 0.31 = 0.0.
pub fn blast_specificity_score(hits: u32, second_hit_pct: Option<f64>) -> f64 {
    // 0 hits: BLAST found no alignments at all — the bait is completely unique;
    // no identity penalty applies.
    if hits == 0 {
        return 1.0;
    }

    let hit_factor = match hits {
        // 1 hit: only the expected self-alignment — also perfectly specific.
        1 => 1.0,
        _ => (1.0 - (hits as f64).log10() / 2.0).clamp(0.0, 1.0),
    };

    let identity_factor = match second_hit_pct {
        Some(pct) => 1.0 - ((pct - 75.0) / 25.0).clamp(0.0, 1.0),
        None => 1.0,
    };

    hit_factor * identity_factor
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::BaitMetric;

    #[test]
    fn test_gc_score_optimal_range() {
        assert_eq!(gc_score(0.40), 1.0);
        assert_eq!(gc_score(0.50), 1.0);
        assert_eq!(gc_score(0.60), 1.0);
    }

    #[test]
    fn test_gc_score_marginal_range() {
        assert_eq!(gc_score(0.35), 0.5);
        assert_eq!(gc_score(0.65), 0.5);
    }

    #[test]
    fn test_gc_score_extreme_range() {
        assert_eq!(gc_score(0.20), 0.0);
        assert_eq!(gc_score(0.80), 0.0);
    }

    #[test]
    fn test_homopolymer_score_no_repeat() {
        // run of 1 → 1.0
        assert_eq!(homopolymer_score(1), 1.0);
    }

    #[test]
    fn test_homopolymer_score_short_run() {
        // run of 2 → still 1.0 (below the 3-base penalty threshold)
        assert_eq!(homopolymer_score(2), 1.0);
    }

    #[test]
    fn test_homopolymer_score_run_of_6() {
        // (6 - 2) / 8 = 0.5 → score = 0.5
        assert!((homopolymer_score(6) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_homopolymer_score_run_of_10_or_more_clamps_to_zero() {
        assert_eq!(homopolymer_score(10), 0.0);
        assert_eq!(homopolymer_score(20), 0.0);
    }

    #[test]
    fn test_mfe_score_zero_is_perfect() {
        assert_eq!(mfe_score(0.0), 1.0);
    }

    #[test]
    fn test_mfe_score_positive_clamps_to_one() {
        assert_eq!(mfe_score(5.0), 1.0);
    }

    #[test]
    fn test_mfe_score_threshold_is_perfect() {
        // Anything at or above -9 kcal/mol scores 1.0 — no penalty.
        assert_eq!(mfe_score(-9.0), 1.0);
        assert_eq!(mfe_score(-5.0), 1.0);
    }

    #[test]
    fn test_mfe_score_at_negative_fifteen_is_zero() {
        assert_eq!(mfe_score(-15.0), 0.0);
    }

    #[test]
    fn test_mfe_score_below_negative_fifteen_clamps_to_zero() {
        assert_eq!(mfe_score(-20.0), 0.0);
        assert_eq!(mfe_score(-50.0), 0.0);
    }

    #[test]
    fn test_mfe_score_midpoint_taper() {
        // At -12 (midpoint of [-9, -15]): t = 0.5, score = (1-0.5)² = 0.25.
        assert!((mfe_score(-12.0) - 0.25).abs() < 1e-9);
    }

    #[test]
    fn test_blast_specificity_zero_hits() {
        assert_eq!(blast_specificity_score(0, None), 1.0);
    }

    #[test]
    fn test_blast_specificity_zero_hits_is_one() {
        // 0 hits means BLAST found no alignments — the bait is completely unique.
        assert_eq!(blast_specificity_score(0, None), 1.0);
        assert_eq!(blast_specificity_score(0, Some(99.0)), 1.0);
    }

    #[test]
    fn test_blast_specificity_one_hit() {
        assert_eq!(blast_specificity_score(1, None), 1.0);
    }

    #[test]
    fn test_blast_specificity_10_hits_is_half() {
        // log10(10) = 1.0 → 1 - 1/2 = 0.5
        assert!((blast_specificity_score(10, None) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_blast_specificity_100_hits_is_zero() {
        // log10(100) = 2.0 → 1 - 2/2 = 0.0
        assert_eq!(blast_specificity_score(100, None), 0.0);
    }

    #[test]
    fn test_blast_specificity_many_hits_clamps_to_zero() {
        assert_eq!(blast_specificity_score(74494, None), 0.0);
    }

    #[test]
    fn test_blast_specificity_identity_no_penalty_below_75_pct() {
        // second hit at 70% → identity_factor = 1.0 (no penalty)
        assert_eq!(blast_specificity_score(1, Some(70.0)), 1.0);
    }

    #[test]
    fn test_blast_specificity_identity_full_penalty_at_100_pct() {
        // second hit at 100% → identity_factor = 0.0
        assert_eq!(blast_specificity_score(1, Some(100.0)), 0.0);
    }

    #[test]
    fn test_blast_specificity_identity_half_penalty_at_87_5_pct() {
        // (87.5 - 75) / 25 = 0.5 → identity_factor = 0.5
        assert!((blast_specificity_score(1, Some(87.5)) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_blast_specificity_rad51_like_bait_is_zero() {
        // 74494 hits → hit_factor = 0.0; identity penalty doesn't matter
        assert_eq!(blast_specificity_score(74494, Some(92.437)), 0.0);
    }

    #[test]
    fn test_bait_score_none_when_no_components() {
        let m = BaitMetric::from_position("b".to_string(), "chr1", 0, 100);
        assert!(bait_score(&m).is_none());
    }

    #[test]
    fn test_bait_score_none_when_one_quality_component_no_blast() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 100);
        m.gc_content = Some(0.50);
        assert!(bait_score(&m).is_none());
    }

    #[test]
    fn test_bait_score_some_when_one_quality_component_plus_blast() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 100);
        m.gc_content = Some(0.50);
        m.blast_hits = Some(1);
        assert!(bait_score(&m).is_some());
    }

    #[test]
    fn test_bait_score_perfect_bait_is_one() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 100);
        m.gc_content = Some(0.50);
        m.longest_homopolymer_size = Some(1);
        m.min_free_energy = Some(0.0);
        m.mean_mappability = Some(1.0);
        m.blast_hits = Some(1);
        let s = bait_score(&m).unwrap();
        assert!((s - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_bait_score_non_specific_bait_is_near_zero() {
        // Good synthesis metrics but 74k BLAST hits → score should be ~0
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 120);
        m.gc_content = Some(0.558);
        m.longest_homopolymer_size = Some(3);
        m.min_free_energy = Some(-15.93);
        m.blast_hits = Some(74494);
        m.blast_second_hit_percent = Some(92.437);
        let s = bait_score(&m).unwrap();
        assert!(s < 0.01, "expected near-zero score, got {s}");
    }

    #[test]
    fn test_bait_score_is_between_zero_and_one() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 120);
        m.gc_content = Some(0.25);
        m.longest_homopolymer_size = Some(8);
        m.blast_hits = Some(5);
        let s = bait_score(&m).unwrap();
        assert!(s >= 0.0);
        assert!(s <= 1.0);
    }

    #[test]
    fn test_bait_score_bait_on_mito_no_mito_identity_penalty() {
        let mut m = BaitMetric::from_position("b".to_string(), "chrM", 0, 120);
        m.gc_content = Some(0.50);
        m.longest_homopolymer_size = Some(2);
        m.mito_hit_identity = Some(100.0);
        let s = bait_score(&m).unwrap();
        assert!(
            s > 0.9,
            "bait on chrM should not be penalized for mito identity, got {s}"
        );
    }

    #[test]
    fn test_bait_score_target_on_mito_no_mito_identity_penalty() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 120);
        m.gc_content = Some(0.50);
        m.longest_homopolymer_size = Some(2);
        m.target_interval = Some("chrM:1-1000".to_string());
        m.mito_hit_identity = Some(100.0);
        let s = bait_score(&m).unwrap();
        assert!(
            s > 0.9,
            "bait targeting chrM should not be penalized for mito identity, got {s}"
        );
    }

    #[test]
    fn test_bait_score_mito_hit_at_100pct_identity_drives_to_zero() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 120);
        m.gc_content = Some(0.50);
        m.longest_homopolymer_size = Some(2);
        m.mito_hit_identity = Some(100.0);
        let s = bait_score(&m).unwrap();
        assert!(
            (s - 0.0).abs() < 1e-9,
            "100% mito identity should drive score to 0, got {s}"
        );
    }

    #[test]
    fn test_bait_score_mito_hit_at_75pct_identity_no_penalty() {
        let mut m_with = BaitMetric::from_position("b".to_string(), "chr1", 0, 120);
        m_with.gc_content = Some(0.50);
        m_with.longest_homopolymer_size = Some(2);
        m_with.mito_hit_identity = Some(75.0);

        let mut m_without = BaitMetric::from_position("b".to_string(), "chr1", 0, 120);
        m_without.gc_content = Some(0.50);
        m_without.longest_homopolymer_size = Some(2);

        let s_with = bait_score(&m_with).unwrap();
        let s_without = bait_score(&m_without).unwrap();
        assert!(
            (s_with - s_without).abs() < 1e-9,
            "75% mito identity should have no penalty"
        );
    }

    #[test]
    fn test_bait_score_mito_hit_partial_penalty() {
        let mut m = BaitMetric::from_position("b".to_string(), "chr1", 0, 120);
        m.gc_content = Some(0.50);
        m.longest_homopolymer_size = Some(2);
        m.mito_hit_identity = Some(87.5);
        let s = bait_score(&m).unwrap();
        let mut m_no_mito = BaitMetric::from_position("b".to_string(), "chr1", 0, 120);
        m_no_mito.gc_content = Some(0.50);
        m_no_mito.longest_homopolymer_size = Some(2);
        let s_no_mito = bait_score(&m_no_mito).unwrap();
        assert!(
            (s - s_no_mito * 0.5).abs() < 1e-9,
            "87.5% identity should apply 50% mito penalty"
        );
    }
}
