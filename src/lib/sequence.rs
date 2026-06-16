//! Sequence-based bait metrics computed directly from the nucleotide string.

use anyhow::{Context, Result};

/// Return the reverse complement of a DNA sequence (case-preserving, IUPAC-aware).
///
/// Delegates to [`bio::alphabets::dna::revcomp`]. All standard IUPAC ambiguity codes
/// (R, Y, S, W, K, M, B, D, H, V, N) and their lowercase counterparts are complemented
/// correctly with casing preserved. Unrecognized ASCII characters are passed through
/// unchanged.
///
/// `bio::alphabets::dna::revcomp` works on raw bytes and reverses them, so a multi-byte
/// (non-ASCII) input would reverse into an invalid UTF-8 byte sequence. Rather than
/// panic, this returns an `Err` in that case. DNA sequences are ASCII, so it never
/// fails for them.
///
/// # Examples
/// ```
/// use chumlib::sequence::reverse_complement;
/// assert_eq!(reverse_complement("ACGT").unwrap(), "ACGT");
/// assert_eq!(reverse_complement("AAACCC").unwrap(), "GGGTTT");
/// assert_eq!(reverse_complement("acgt").unwrap(), "acgt");
/// ```
pub fn reverse_complement(seq: &str) -> Result<String> {
    String::from_utf8(bio::alphabets::dna::revcomp(seq.as_bytes())).with_context(
        || "reverse complement produced non-UTF-8 output; the sequence contains non-ASCII bytes",
    )
}

/// Compute GC content as a fraction.
///
/// `N` and `.` (no-call / gap) bases are excluded from the denominator. Lowercase
/// soft-masked bases are case-folded and counted like their uppercase forms; they are
/// **not** excluded (despite being "masked" in the sense of [`masked_bases`]).
///
/// # Examples
/// ```
/// use chumlib::sequence::gc_content;
/// assert!((gc_content("ACGT") - 0.5).abs() < 1e-9);
/// assert!((gc_content("NNNN") - 0.0).abs() < 1e-9);
/// ```
pub fn gc_content(seq: &str) -> f64 {
    let countable: Vec<char> = seq
        .chars()
        .filter(|c| !matches!(c.to_ascii_uppercase(), 'N' | '.'))
        .collect();
    let total = countable.len();
    if total == 0 {
        return 0.0;
    }
    let gc = countable
        .iter()
        .filter(|c| matches!(c.to_ascii_uppercase(), 'G' | 'C'))
        .count();
    gc as f64 / total as f64
}

/// Count masked bases: lowercase letters, `N`, and `.`.
///
/// # Examples
/// ```
/// use chumlib::sequence::masked_bases;
/// assert_eq!(masked_bases("ACGTacgtNN.."), 8);
/// assert_eq!(masked_bases("ACGT"), 0);
/// ```
pub fn masked_bases(seq: &str) -> u32 {
    seq.chars()
        .filter(|c| c.is_lowercase() || *c == 'N' || *c == '.')
        .count() as u32
}

/// Count the number of maximal homopolymer runs of length ≥ `min_len` (case-insensitive).
///
/// Only runs of the four nucleotides A/C/G/T are counted. Runs of `N`, `.`, gaps, or
/// any other character are not homopolymers and are ignored (they still break a run).
///
/// # Examples
/// ```
/// use chumlib::sequence::count_homopolymers_min;
/// // ACGTTTTCAAAA has runs: TTTT (len 4) and AAAA (len 4)
/// assert_eq!(count_homopolymers_min("ACGTTTTCAAAA", 3), 2);
/// assert_eq!(count_homopolymers_min("ACGT", 3), 0);
/// assert_eq!(count_homopolymers_min("AAACCC", 3), 2);
/// // A run of N is not a homopolymer.
/// assert_eq!(count_homopolymers_min("ACNNNNNGT", 3), 0);
/// ```
pub fn count_homopolymers_min(seq: &str, min_len: usize) -> u32 {
    if seq.is_empty() {
        return 0;
    }
    let chars: Vec<char> = seq.chars().map(|c| c.to_ascii_uppercase()).collect();
    let mut count = 0u32;
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        let mut j = i + 1;
        while j < chars.len() && chars[j] == c {
            j += 1;
        }
        if matches!(c, 'A' | 'C' | 'G' | 'T') && j - i >= min_len {
            count += 1;
        }
        i = j;
    }
    count
}

/// Count distinct homopolymer runs of length ≥ 3 (case-insensitive).
///
/// # Examples
/// ```
/// use chumlib::sequence::homopolymers_size_3_or_greater;
/// assert_eq!(homopolymers_size_3_or_greater("ACGTTTTCAAAA"), 2);
/// assert_eq!(homopolymers_size_3_or_greater("ACGT"), 0);
/// ```
pub fn homopolymers_size_3_or_greater(seq: &str) -> u32 {
    count_homopolymers_min(seq, 3)
}

/// Return `true` when `name` looks like a mitochondrial sequence identifier.
///
/// Recognizes the common names used across reference builds:
/// `chrM`, `MT`, `M`, `chrMT`, `NC_012920` (and versioned suffixes such as
/// `NC_012920.1`), and any name containing the substring `mitochon`.
pub fn is_mito_chrom(name: &str) -> bool {
    matches!(name, "chrM" | "MT" | "M" | "chrMT")
        || name.starts_with("NC_012920")
        || name.to_ascii_lowercase().contains("mitochon")
}

/// Return the length of the longest homopolymer run of A/C/G/T (case-insensitive).
///
/// Runs of `N`, `.`, gaps, or any non-nucleotide character are not homopolymers and are
/// ignored. Returns 0 for empty sequences (or sequences with no A/C/G/T) and 1 when no
/// two consecutive nucleotides are identical.
///
/// # Examples
/// ```
/// use chumlib::sequence::longest_homopolymer_size;
/// assert_eq!(longest_homopolymer_size("ACGTTTTCAAAA"), 4);
/// assert_eq!(longest_homopolymer_size("ACGT"), 1);
/// assert_eq!(longest_homopolymer_size(""), 0);
/// // A run of N does not count as a homopolymer.
/// assert_eq!(longest_homopolymer_size("ACNNNNNGT"), 1);
/// ```
pub fn longest_homopolymer_size(seq: &str) -> u32 {
    if seq.is_empty() {
        return 0;
    }
    let chars: Vec<char> = seq.chars().map(|c| c.to_ascii_uppercase()).collect();
    let mut max_run = 0u32;
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        let mut j = i + 1;
        while j < chars.len() && chars[j] == c {
            j += 1;
        }
        let run_len = (j - i) as u32;
        if matches!(c, 'A' | 'C' | 'G' | 'T') && run_len > max_run {
            max_run = run_len;
        }
        i = j;
    }
    max_run
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_mito_chrom_various_names() {
        assert!(is_mito_chrom("chrM"));
        assert!(is_mito_chrom("MT"));
        assert!(is_mito_chrom("M"));
        assert!(is_mito_chrom("chrMT"));
        assert!(is_mito_chrom("NC_012920"));
        assert!(is_mito_chrom("NC_012920.1"));
        assert!(is_mito_chrom("human_mitochondrion"));
        assert!(!is_mito_chrom("chr1"));
        assert!(!is_mito_chrom("chrX"));
        assert!(!is_mito_chrom("MT_random"));
    }

    #[test]
    fn test_gc_content_even_split() {
        // ACGT: G+C = 2 out of 4 countable
        assert!((gc_content("ACGT") - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_gc_content_all_gc() {
        assert!((gc_content("GCGC") - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_gc_content_all_at() {
        assert_eq!(gc_content("ATAT"), 0.0);
    }

    #[test]
    fn test_gc_content_excludes_n_from_denominator() {
        // "ACGN": countable = ACG (3 bases), G+C = 2 → 2/3
        assert!((gc_content("ACGN") - 2.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_gc_content_all_n_is_zero() {
        assert_eq!(gc_content("NNNN"), 0.0);
    }

    #[test]
    fn test_gc_content_empty_is_zero() {
        assert_eq!(gc_content(""), 0.0);
    }

    #[test]
    fn test_gc_content_lowercase_counts_as_gc() {
        // lowercase gc should still count as G+C
        assert!((gc_content("acgt") - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_gc_content_known_bait() {
        // "ACGTTTTCAAAA" → GC = 3/12 = 0.25
        assert!((gc_content("ACGTTTTCAAAA") - 0.25).abs() < 1e-9);
    }

    #[test]
    fn test_gc_content_gc_rich_bait() {
        // "ACAGGGG" → GC = 5/7
        assert!((gc_content("ACAGGGG") - 5.0 / 7.0).abs() < 1e-9);
    }

    #[test]
    fn test_masked_bases_none() {
        assert_eq!(masked_bases("ACGT"), 0);
    }

    #[test]
    fn test_masked_bases_lowercase() {
        assert_eq!(masked_bases("acgt"), 4);
    }

    #[test]
    fn test_masked_bases_n_and_dot() {
        assert_eq!(masked_bases("NN.."), 4);
    }

    #[test]
    fn test_masked_bases_mixed() {
        assert_eq!(masked_bases("ACGTacgtNN.."), 8);
    }

    #[test]
    fn test_masked_bases_empty() {
        assert_eq!(masked_bases(""), 0);
    }

    #[test]
    fn test_count_homopolymers_min_none() {
        assert_eq!(count_homopolymers_min("ACGT", 3), 0);
    }

    #[test]
    fn test_count_homopolymers_min_one_run() {
        assert_eq!(count_homopolymers_min("ACGTTTACGT", 3), 1);
    }

    #[test]
    fn test_count_homopolymers_min_two_runs() {
        // TTTT and AAAA → 2
        assert_eq!(count_homopolymers_min("ACGTTTTCAAAA", 3), 2);
    }

    #[test]
    fn test_count_homopolymers_min_run_exactly_at_boundary() {
        assert_eq!(count_homopolymers_min("AAA", 3), 1);
        assert_eq!(count_homopolymers_min("AA", 3), 0);
    }

    #[test]
    fn test_count_homopolymers_min_case_insensitive() {
        // "acgtttacgt" = lowercase; "ttt" → 1 run
        assert_eq!(count_homopolymers_min("acgtttacgt", 3), 1);
    }

    #[test]
    fn test_count_homopolymers_min_empty() {
        assert_eq!(count_homopolymers_min("", 3), 0);
    }

    #[test]
    fn test_count_homopolymers_min_acagggg() {
        // "ACAGGGG": GGGG is a run of 4 → 1 run ≥ 3
        assert_eq!(count_homopolymers_min("ACAGGGG", 3), 1);
    }

    #[test]
    fn test_homopolymers_size_3_or_greater_two_runs() {
        // "ACGTTTTCAAAA": TTTT (4) and AAAA (4) → 2
        assert_eq!(homopolymers_size_3_or_greater("ACGTTTTCAAAA"), 2);
    }

    #[test]
    fn test_homopolymers_size_3_or_greater_none() {
        assert_eq!(homopolymers_size_3_or_greater("ACGT"), 0);
    }

    #[test]
    fn test_longest_homopolymer_empty() {
        assert_eq!(longest_homopolymer_size(""), 0);
    }

    #[test]
    fn test_longest_homopolymer_no_repeat() {
        assert_eq!(longest_homopolymer_size("ACGT"), 1);
    }

    #[test]
    fn test_longest_homopolymer_two_equal_runs() {
        // "ACGTTTTCAAAA": TTTT=4, AAAA=4 → max = 4
        assert_eq!(longest_homopolymer_size("ACGTTTTCAAAA"), 4);
    }

    #[test]
    fn test_longest_homopolymer_gc_rich_bait() {
        // "ACAGGGG": GGGG = 4
        assert_eq!(longest_homopolymer_size("ACAGGGG"), 4);
    }

    #[test]
    fn test_reverse_complement_palindrome() {
        assert_eq!(reverse_complement("ACGT").unwrap(), "ACGT");
    }

    #[test]
    fn test_reverse_complement_simple() {
        assert_eq!(reverse_complement("AAACCC").unwrap(), "GGGTTT");
    }

    #[test]
    fn test_reverse_complement_preserves_case() {
        // Lowercase (repeat-masked) bases stay lowercase after complementing.
        assert_eq!(reverse_complement("acgt").unwrap(), "acgt");
        // Mixed case: reversed "AcGt" → t,G,c,A → complement a,C,g,T → "aCgT"
        assert_eq!(reverse_complement("AcGt").unwrap(), "aCgT");
    }

    #[test]
    fn test_reverse_complement_n_passthrough() {
        assert_eq!(reverse_complement("ACNGT").unwrap(), "ACNGT");
    }

    #[test]
    fn test_reverse_complement_homopolymer_conversion() {
        // ACAGGGG on + strand → probe on - strand is RC: CCCCTGT
        assert_eq!(reverse_complement("ACAGGGG").unwrap(), "CCCCTGT");
    }

    #[test]
    fn test_longest_homopolymer_reverse_complement() {
        // "TCCCCTG": CCCC = 4
        assert_eq!(longest_homopolymer_size("TCCCCTG"), 4);
    }

    #[test]
    fn test_longest_homopolymer_case_insensitive() {
        assert_eq!(longest_homopolymer_size("acgttttcaaaa"), 4);
    }

    #[test]
    fn test_longest_homopolymer_whole_sequence_same() {
        assert_eq!(longest_homopolymer_size("AAAAAAA"), 7);
    }

    #[test]
    fn test_reverse_complement_non_ascii_returns_err_not_panic() {
        // Multi-byte UTF-8 reverses into invalid UTF-8; must return Err, not panic.
        assert!(reverse_complement("é").is_err());
        // The U+FFFD replacement char (what from_utf8_lossy emits for bad bytes).
        assert!(reverse_complement("AC\u{FFFD}GT").is_err());
    }

    #[test]
    fn test_longest_homopolymer_ignores_n_run() {
        // A run of N is an assembly gap, not a homopolymer.
        assert_eq!(longest_homopolymer_size("ACNNNNNGT"), 1);
        assert_eq!(longest_homopolymer_size("NNNNN"), 0);
        // A real homopolymer adjacent to an N run is still counted.
        assert_eq!(longest_homopolymer_size("AAAANNNNN"), 4);
    }

    #[test]
    fn test_count_homopolymers_ignores_n_run() {
        assert_eq!(count_homopolymers_min("ACNNNNNGT", 3), 0);
        // Real run plus an N run → only the real run counts.
        assert_eq!(count_homopolymers_min("AAAANNNNN", 3), 1);
    }
}
