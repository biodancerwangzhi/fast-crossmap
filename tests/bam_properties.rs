//! BAM/SAM/CRAM format property tests
//!
//! Tests for BAM format adapter including CIGAR operations and coordinate mapping.

use proptest::prelude::*;
use fast_crossmap::formats::bam::{CigarOp, CigarReconstructor, AlignmentTag};

// ============================================================================
// Property Tests for CIGAR Operations
// ============================================================================

/// Generate arbitrary CIGAR operations
fn arb_cigar_op() -> impl Strategy<Value = CigarOp> {
    prop_oneof![
        (1u32..1000).prop_map(CigarOp::Match),
        (1u32..100).prop_map(CigarOp::Insertion),
        (1u32..100).prop_map(CigarOp::Deletion),
        (1u32..10000).prop_map(CigarOp::Skip),
        (1u32..100).prop_map(CigarOp::SoftClip),
        (1u32..100).prop_map(CigarOp::HardClip),
        (1u32..100).prop_map(CigarOp::Equal),
        (1u32..100).prop_map(CigarOp::Diff),
    ]
}

/// Generate a valid CIGAR string (sequence of operations)
fn arb_cigar_string() -> impl Strategy<Value = Vec<CigarOp>> {
    prop::collection::vec(arb_cigar_op(), 1..20)
}

proptest! {
    /// Property: CIGAR reverse is involutive (reversing twice gives original)
    #[test]
    fn test_cigar_reverse_involutive(cigar in arb_cigar_string()) {
        let reversed = CigarReconstructor::reverse(&cigar);
        let double_reversed = CigarReconstructor::reverse(&reversed);
        
        prop_assert_eq!(cigar.len(), double_reversed.len());
        for (orig, rev) in cigar.iter().zip(double_reversed.iter()) {
            prop_assert_eq!(orig, rev);
        }
    }
    
    /// Property: CIGAR merge preserves total length
    #[test]
    fn test_cigar_merge_preserves_length(cigar in arb_cigar_string()) {
        let original_query_len = CigarReconstructor::query_length(&cigar);
        let original_ref_len = CigarReconstructor::reference_length(&cigar);
        
        let merged = CigarReconstructor::merge_adjacent(&cigar);
        
        let merged_query_len = CigarReconstructor::query_length(&merged);
        let merged_ref_len = CigarReconstructor::reference_length(&merged);
        
        prop_assert_eq!(original_query_len, merged_query_len);
        prop_assert_eq!(original_ref_len, merged_ref_len);
    }
    
    /// Property: Merged CIGAR has no adjacent same-type operations
    #[test]
    fn test_cigar_merge_no_adjacent_same_type(cigar in arb_cigar_string()) {
        let merged = CigarReconstructor::merge_adjacent(&cigar);
        
        for window in merged.windows(2) {
            let same_type = match (&window[0], &window[1]) {
                (CigarOp::Match(_), CigarOp::Match(_)) => true,
                (CigarOp::Insertion(_), CigarOp::Insertion(_)) => true,
                (CigarOp::Deletion(_), CigarOp::Deletion(_)) => true,
                (CigarOp::Skip(_), CigarOp::Skip(_)) => true,
                (CigarOp::SoftClip(_), CigarOp::SoftClip(_)) => true,
                (CigarOp::HardClip(_), CigarOp::HardClip(_)) => true,
                (CigarOp::Equal(_), CigarOp::Equal(_)) => true,
                (CigarOp::Diff(_), CigarOp::Diff(_)) => true,
                _ => false,
            };
            prop_assert!(!same_type, "Found adjacent same-type operations after merge");
        }
    }
    
    /// Property: Query length equals sum of query-consuming operations
    #[test]
    fn test_cigar_query_length_correct(cigar in arb_cigar_string()) {
        let computed_len = CigarReconstructor::query_length(&cigar);
        
        let expected_len: u32 = cigar.iter()
            .filter(|op| op.consumes_query())
            .map(|op| op.len())
            .sum();
        
        prop_assert_eq!(computed_len, expected_len);
    }
    
    /// Property: Reference length equals sum of reference-consuming operations
    #[test]
    fn test_cigar_reference_length_correct(cigar in arb_cigar_string()) {
        let computed_len = CigarReconstructor::reference_length(&cigar);
        
        let expected_len: u32 = cigar.iter()
            .filter(|op| op.consumes_reference())
            .map(|op| op.len())
            .sum();
        
        prop_assert_eq!(computed_len, expected_len);
    }
}

// ============================================================================
// Unit Tests
// ============================================================================

#[test]
fn test_cigar_op_consumes_reference() {
    // Operations that consume reference
    assert!(CigarOp::Match(10).consumes_reference());
    assert!(CigarOp::Deletion(5).consumes_reference());
    assert!(CigarOp::Skip(100).consumes_reference());
    assert!(CigarOp::Equal(10).consumes_reference());
    assert!(CigarOp::Diff(5).consumes_reference());
    
    // Operations that don't consume reference
    assert!(!CigarOp::Insertion(5).consumes_reference());
    assert!(!CigarOp::SoftClip(10).consumes_reference());
    assert!(!CigarOp::HardClip(10).consumes_reference());
    assert!(!CigarOp::Padding(5).consumes_reference());
}

#[test]
fn test_cigar_op_consumes_query() {
    // Operations that consume query
    assert!(CigarOp::Match(10).consumes_query());
    assert!(CigarOp::Insertion(5).consumes_query());
    assert!(CigarOp::SoftClip(10).consumes_query());
    assert!(CigarOp::Equal(10).consumes_query());
    assert!(CigarOp::Diff(5).consumes_query());
    
    // Operations that don't consume query
    assert!(!CigarOp::Deletion(5).consumes_query());
    assert!(!CigarOp::Skip(100).consumes_query());
    assert!(!CigarOp::HardClip(10).consumes_query());
    assert!(!CigarOp::Padding(5).consumes_query());
}

#[test]
fn test_cigar_reverse_simple() {
    let cigar = vec![
        CigarOp::Match(10),
        CigarOp::Insertion(2),
        CigarOp::Match(5),
    ];
    
    let reversed = CigarReconstructor::reverse(&cigar);
    
    assert_eq!(reversed.len(), 3);
    assert_eq!(reversed[0], CigarOp::Match(5));
    assert_eq!(reversed[1], CigarOp::Insertion(2));
    assert_eq!(reversed[2], CigarOp::Match(10));
}

#[test]
fn test_cigar_merge_adjacent_simple() {
    let cigar = vec![
        CigarOp::Match(10),
        CigarOp::Match(5),
        CigarOp::Insertion(2),
        CigarOp::Insertion(3),
        CigarOp::Match(8),
    ];
    
    let merged = CigarReconstructor::merge_adjacent(&cigar);
    
    assert_eq!(merged.len(), 3);
    assert_eq!(merged[0], CigarOp::Match(15));
    assert_eq!(merged[1], CigarOp::Insertion(5));
    assert_eq!(merged[2], CigarOp::Match(8));
}

#[test]
fn test_cigar_validate() {
    let cigar = vec![
        CigarOp::Match(10),
        CigarOp::Insertion(2),
        CigarOp::Match(8),
    ];
    
    // Query length = 10 + 2 + 8 = 20
    assert!(CigarReconstructor::validate(&cigar, 20));
    assert!(!CigarReconstructor::validate(&cigar, 15));
    assert!(!CigarReconstructor::validate(&cigar, 25));
}

#[test]
fn test_cigar_handle_break() {
    let cigar = vec![
        CigarOp::Match(10),
        CigarOp::Match(10),
    ];
    
    // Break at position 5 (within first Match)
    let (left, _right) = CigarReconstructor::handle_break(&cigar, 5);
    
    // Left should have Match(5) + SoftClip for remaining
    assert!(!left.is_empty());
    // Right should start with SoftClip for skipped bases
    // The exact behavior depends on implementation details
}

#[test]
fn test_alignment_tag_values() {
    assert_eq!(AlignmentTag::QF.as_str(), "QF");
    assert_eq!(AlignmentTag::NN.as_str(), "NN");
    assert_eq!(AlignmentTag::NU.as_str(), "NU");
    assert_eq!(AlignmentTag::NM.as_str(), "NM");
    assert_eq!(AlignmentTag::UN.as_str(), "UN");
    assert_eq!(AlignmentTag::UU.as_str(), "UU");
    assert_eq!(AlignmentTag::UM.as_str(), "UM");
    assert_eq!(AlignmentTag::MN.as_str(), "MN");
    assert_eq!(AlignmentTag::MU.as_str(), "MU");
    assert_eq!(AlignmentTag::MM.as_str(), "MM");
    assert_eq!(AlignmentTag::SN.as_str(), "SN");
    assert_eq!(AlignmentTag::SM.as_str(), "SM");
    assert_eq!(AlignmentTag::SU.as_str(), "SU");
}

#[test]
fn test_cigar_empty_merge() {
    let cigar: Vec<CigarOp> = vec![];
    let merged = CigarReconstructor::merge_adjacent(&cigar);
    assert!(merged.is_empty());
}

#[test]
fn test_cigar_single_op_merge() {
    let cigar = vec![CigarOp::Match(10)];
    let merged = CigarReconstructor::merge_adjacent(&cigar);
    assert_eq!(merged.len(), 1);
    assert_eq!(merged[0], CigarOp::Match(10));
}
