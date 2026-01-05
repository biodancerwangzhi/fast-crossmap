//! DNA sequence utilities
//!
//! Provides functions for DNA sequence manipulation including
//! reverse complement and validation.

/// Complement a single DNA base
/// 
/// Supports standard bases (A, T, G, C) and IUPAC ambiguity codes.
/// Returns the same character for non-DNA characters.
#[inline]
pub fn complement_base(base: u8) -> u8 {
    match base {
        // Standard bases
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'a' => b't',
        b't' => b'a',
        b'g' => b'c',
        b'c' => b'g',
        // IUPAC ambiguity codes
        b'R' => b'Y', // R = A or G -> Y = T or C
        b'Y' => b'R', // Y = T or C -> R = A or G
        b'S' => b'S', // S = G or C -> S = G or C (self-complementary)
        b'W' => b'W', // W = A or T -> W = A or T (self-complementary)
        b'K' => b'M', // K = G or T -> M = A or C
        b'M' => b'K', // M = A or C -> K = G or T
        b'B' => b'V', // B = C, G, or T -> V = A, C, or G
        b'V' => b'B', // V = A, C, or G -> B = C, G, or T
        b'D' => b'H', // D = A, G, or T -> H = A, C, or T
        b'H' => b'D', // H = A, C, or T -> D = A, G, or T
        b'N' => b'N', // N = any base -> N = any base
        b'r' => b'y',
        b'y' => b'r',
        b's' => b's',
        b'w' => b'w',
        b'k' => b'm',
        b'm' => b'k',
        b'b' => b'v',
        b'v' => b'b',
        b'd' => b'h',
        b'h' => b'd',
        b'n' => b'n',
        // Unknown - return as-is
        _ => base,
    }
}

/// Compute the reverse complement of a DNA sequence
/// 
/// # Arguments
/// * `seq` - DNA sequence as a string slice
/// 
/// # Returns
/// Reverse complement as a new String
/// 
/// # Examples
/// ```
/// use fast_crossmap::core::dna::revcomp;
/// 
/// assert_eq!(revcomp("AACGT"), "ACGTT");
/// assert_eq!(revcomp("ATGC"), "GCAT");
/// assert_eq!(revcomp(""), "");
/// ```
pub fn revcomp(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(complement_base)
        .map(|b| b as char)
        .collect()
}

/// Compute the reverse complement of a DNA sequence in-place (bytes)
/// 
/// # Arguments
/// * `seq` - Mutable byte slice to reverse complement in place
pub fn revcomp_inplace(seq: &mut [u8]) {
    seq.reverse();
    for base in seq.iter_mut() {
        *base = complement_base(*base);
    }
}

/// Check if a character is a valid DNA base (standard or IUPAC)
/// 
/// # Arguments
/// * `base` - Character to check
/// 
/// # Returns
/// true if the character is a valid DNA base
#[inline]
pub fn is_dna_base(base: u8) -> bool {
    matches!(
        base,
        b'A' | b'T' | b'G' | b'C' | b'a' | b't' | b'g' | b'c' |
        b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'V' | b'D' | b'H' | b'N' |
        b'r' | b'y' | b's' | b'w' | b'k' | b'm' | b'b' | b'v' | b'd' | b'h' | b'n'
    )
}

/// Check if a string is a valid DNA sequence
/// 
/// # Arguments
/// * `seq` - String to check
/// 
/// # Returns
/// true if all characters are valid DNA bases
/// 
/// # Examples
/// ```
/// use fast_crossmap::core::dna::is_dna;
/// 
/// assert!(is_dna("ATGC"));
/// assert!(is_dna("atgc"));
/// assert!(is_dna("ATGCN")); // N is valid IUPAC
/// assert!(!is_dna("ATGCX")); // X is not valid
/// assert!(is_dna("")); // Empty is valid
/// ```
pub fn is_dna(seq: &str) -> bool {
    seq.bytes().all(is_dna_base)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement_base_standard() {
        assert_eq!(complement_base(b'A'), b'T');
        assert_eq!(complement_base(b'T'), b'A');
        assert_eq!(complement_base(b'G'), b'C');
        assert_eq!(complement_base(b'C'), b'G');
    }

    #[test]
    fn test_complement_base_lowercase() {
        assert_eq!(complement_base(b'a'), b't');
        assert_eq!(complement_base(b't'), b'a');
        assert_eq!(complement_base(b'g'), b'c');
        assert_eq!(complement_base(b'c'), b'g');
    }

    #[test]
    fn test_complement_base_iupac() {
        assert_eq!(complement_base(b'R'), b'Y');
        assert_eq!(complement_base(b'Y'), b'R');
        assert_eq!(complement_base(b'N'), b'N');
    }

    #[test]
    fn test_revcomp_basic() {
        assert_eq!(revcomp("AACGT"), "ACGTT");
        assert_eq!(revcomp("ATGC"), "GCAT");
        assert_eq!(revcomp("A"), "T");
        assert_eq!(revcomp(""), "");
    }

    #[test]
    fn test_revcomp_lowercase() {
        assert_eq!(revcomp("aacgt"), "acgtt");
    }

    #[test]
    fn test_revcomp_mixed_case() {
        assert_eq!(revcomp("AaCcGgTt"), "aAcCgGtT");
    }

    #[test]
    fn test_revcomp_inplace() {
        let mut seq = b"ATGC".to_vec();
        revcomp_inplace(&mut seq);
        assert_eq!(seq, b"GCAT");
    }

    #[test]
    fn test_is_dna() {
        assert!(is_dna("ATGC"));
        assert!(is_dna("atgc"));
        assert!(is_dna("ATGCatgc"));
        assert!(is_dna("ATGCN"));
        assert!(is_dna("RYSWKMBVDHN"));
        assert!(is_dna(""));
        assert!(!is_dna("ATGCX"));
        assert!(!is_dna("ATGC "));
        assert!(!is_dna("123"));
    }

    #[test]
    fn test_is_dna_base() {
        assert!(is_dna_base(b'A'));
        assert!(is_dna_base(b'a'));
        assert!(is_dna_base(b'N'));
        assert!(!is_dna_base(b'X'));
        assert!(!is_dna_base(b' '));
    }
}
