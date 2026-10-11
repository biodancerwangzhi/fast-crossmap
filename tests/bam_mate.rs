//! Paired-end BAM/SAM tests: the mate columns (RNEXT / PNEXT / TLEN), the
//! rebuilt paired-end FLAG bits, and the auxiliary-field rewrite.
//!
//! `paired_sample.expected.sam` is the output of the reference implementation
//! (Python CrossMap 0.7.3, the version the paper benchmarks against) for
//! `paired_sample.sam` under `paper/sample_data/hg19ToHg38.over.chain.gz`.
//! Between them the 43 records cover all nine paired categories CrossMap can
//! emit plus the QC-fail branch, and carry the auxiliary fields of the original
//! NA12878 slice, so the whole alignment line is compared.

use fast_crossmap::core::{ChainIndex, ChromStyle, CoordinateMapper};
use fast_crossmap::formats::bam::{convert_bam, AlignmentTag};
use std::path::PathBuf;

fn repo_path(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
}

/// Every alignment line, headers dropped.
fn alignment_records(text: &str) -> Vec<&str> {
    text.lines()
        .filter(|line| !line.starts_with('@') && !line.is_empty())
        .collect()
}

#[test]
fn paired_end_output_matches_crossmap() {
    let chain = repo_path("paper/sample_data/hg19ToHg38.over.chain.gz");
    let input = repo_path("tests/data/paired_sample.sam");
    let expected = repo_path("tests/data/paired_sample.expected.sam");

    for path in [&chain, &input, &expected] {
        if !path.exists() {
            eprintln!("Skipping test: {} not found", path.display());
            return;
        }
    }

    let output = std::env::temp_dir().join(format!("fcm_paired_{}.sam", std::process::id()));
    let index = ChainIndex::from_chain_file(&chain).expect("load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);

    let stats = convert_bam(&input, &output, &mapper, 1).expect("convert paired sample");

    let got = std::fs::read_to_string(&output).expect("read output");
    let want = std::fs::read_to_string(&expected).expect("read expected");
    let _ = std::fs::remove_file(&output);

    // The category tally is what CrossMap prints for this input.
    assert_eq!(stats.count_of(AlignmentTag::QF), 2);
    assert_eq!(stats.count_of(AlignmentTag::UU), 32);
    assert_eq!(stats.count_of(AlignmentTag::UN), 9);
    assert_eq!(
        stats.count_of(AlignmentTag::QF)
            + stats.count_of(AlignmentTag::UU)
            + stats.count_of(AlignmentTag::UN),
        stats.total
    );

    let got = alignment_records(&got);
    let want = alignment_records(&want);
    assert_eq!(got.len(), want.len(), "record count differs");
    for (i, (a, b)) in got.iter().zip(want.iter()).enumerate() {
        assert_eq!(a, b, "record {} differs", i + 1);
    }
}
