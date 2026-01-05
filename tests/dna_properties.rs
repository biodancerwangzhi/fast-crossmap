//! Property-based tests for DNA utilities
//!
//! **Feature: fast-crossmap, Property 8: DNA 反向互补幂等性**
//! **Validates: Requirements 5.2, 7.4, 10.3**

use fast_crossmap::core::dna::{complement_base, is_dna, revcomp, revcomp_inplace};
use proptest::prelude::*;

/// Generate a random DNA sequence using standard bases
fn dna_sequence_strategy() -> impl Strategy<Value = String> {
    prop::collection::vec(prop::sample::select(vec!['A', 'T', 'G', 'C']), 0..100)
        .prop_map(|chars| chars.into_iter().collect())
}

/// Generate a random DNA sequence including IUPAC codes
fn dna_iupac_strategy() -> impl Strategy<Value = String> {
    prop::collection::vec(
        prop::sample::select(vec![
            'A', 'T', 'G', 'C', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'V', 'D', 'H', 'N',
        ]),
        0..100,
    )
    .prop_map(|chars| chars.into_iter().collect())
}

/// Generate a random DNA sequence with mixed case
fn dna_mixed_case_strategy() -> impl Strategy<Value = String> {
    prop::collection::vec(
        prop::sample::select(vec![
            'A', 'T', 'G', 'C', 'a', 't', 'g', 'c',
        ]),
        0..100,
    )
    .prop_map(|chars| chars.into_iter().collect())
}

proptest! {
    /// **Property 8: DNA 反向互补幂等性**
    /// 
    /// For any DNA sequence, applying reverse complement twice should return
    /// the original sequence: revcomp(revcomp(seq)) == seq
    /// 
    /// This is a fundamental property of DNA complementarity.
    #[test]
    fn prop_revcomp_involution(seq in dna_sequence_strategy()) {
        // revcomp(revcomp(x)) == x
        let result = revcomp(&revcomp(&seq));
        prop_assert_eq!(result, seq, "revcomp should be an involution");
    }

    /// Property: Reverse complement with IUPAC codes is also an involution
    #[test]
    fn prop_revcomp_involution_iupac(seq in dna_iupac_strategy()) {
        let result = revcomp(&revcomp(&seq));
        prop_assert_eq!(result, seq, "revcomp with IUPAC codes should be an involution");
    }

    /// Property: Reverse complement preserves case pattern
    #[test]
    fn prop_revcomp_involution_mixed_case(seq in dna_mixed_case_strategy()) {
        let result = revcomp(&revcomp(&seq));
        prop_assert_eq!(result, seq, "revcomp should preserve case");
    }

    /// Property: Reverse complement preserves length
    #[test]
    fn prop_revcomp_preserves_length(seq in dna_sequence_strategy()) {
        let result = revcomp(&seq);
        prop_assert_eq!(result.len(), seq.len(), "revcomp should preserve length");
    }

    /// Property: Complement of complement is identity for single bases
    #[test]
    fn prop_complement_involution(base in prop::sample::select(vec![
        b'A', b'T', b'G', b'C', b'R', b'Y', b'S', b'W', b'K', b'M', b'N'
    ])) {
        let result = complement_base(complement_base(base));
        prop_assert_eq!(result, base, "complement should be an involution");
    }

    /// Property: is_dna returns true for all generated DNA sequences
    #[test]
    fn prop_generated_dna_is_valid(seq in dna_iupac_strategy()) {
        prop_assert!(is_dna(&seq), "Generated DNA sequence should be valid");
    }

    /// Property: revcomp of valid DNA is still valid DNA
    #[test]
    fn prop_revcomp_preserves_validity(seq in dna_iupac_strategy()) {
        let result = revcomp(&seq);
        prop_assert!(is_dna(&result), "revcomp should produce valid DNA");
    }

    /// Property: In-place revcomp produces same result as functional revcomp
    #[test]
    fn prop_revcomp_inplace_equivalent(seq in dna_sequence_strategy()) {
        let functional_result = revcomp(&seq);
        let mut inplace_seq = seq.as_bytes().to_vec();
        revcomp_inplace(&mut inplace_seq);
        let inplace_result = String::from_utf8(inplace_seq).unwrap();
        prop_assert_eq!(functional_result, inplace_result, 
            "In-place and functional revcomp should produce same result");
    }

    /// Property: Empty sequence reverse complement is empty
    #[test]
    fn prop_revcomp_empty(_dummy in Just(())) {
        prop_assert_eq!(revcomp(""), "", "revcomp of empty should be empty");
    }

    /// Property: Single base reverse complement is its complement
    #[test]
    fn prop_revcomp_single_base(base in prop::sample::select(vec!['A', 'T', 'G', 'C'])) {
        let seq = base.to_string();
        let result = revcomp(&seq);
        let expected = match base {
            'A' => "T",
            'T' => "A",
            'G' => "C",
            'C' => "G",
            _ => unreachable!(),
        };
        prop_assert_eq!(result, expected, "Single base revcomp should be complement");
    }
}

/// Additional edge case tests (not property-based)
#[cfg(test)]
mod edge_cases {
    use super::*;

    #[test]
    fn test_revcomp_palindrome() {
        // Some sequences are their own reverse complement (palindromes)
        assert_eq!(revcomp("GCGC"), "GCGC");
        assert_eq!(revcomp("ATAT"), "ATAT");
        assert_eq!(revcomp("CATG"), "CATG");
    }

    #[test]
    fn test_revcomp_known_sequences() {
        // Known biological sequences
        assert_eq!(revcomp("ATG"), "CAT"); // Start codon
        assert_eq!(revcomp("TAA"), "TTA"); // Stop codon
        assert_eq!(revcomp("GAATTC"), "GAATTC"); // EcoRI site (palindrome)
    }
}
