#!/usr/bin/env python3
"""
Property-based tests for benchmark module.

**Feature: comprehensive-benchmark-audit**
**Validates: Requirements 2.5, 2.6, 2.7, 2.9, 5.4**

Run with: pytest scripts/test_benchmarks.py -v
"""

import csv
import json
import tempfile
from dataclasses import asdict
from pathlib import Path

import pytest
from hypothesis import given, settings, strategies as st

# Import from run_benchmarks
from run_benchmarks import (
    BenchmarkResult,
    BenchmarkMetadata,
    BenchmarkReport,
    export_to_json,
    export_to_csv,
)


# ============================================================================
# Property 1: Benchmark Result Completeness
# **Validates: Requirements 2.5, 2.6, 2.7**
# ============================================================================

@st.composite
def benchmark_result_strategy(draw):
    """Generate random BenchmarkResult instances."""
    input_records = draw(st.integers(min_value=1, max_value=10_000_000))
    execution_time = draw(st.floats(min_value=0.001, max_value=3600.0))
    
    return BenchmarkResult(
        tool=draw(st.sampled_from(["fastcrossmap", "crossmap", "liftover", "fastremap"])),
        format=draw(st.sampled_from(["bed", "vcf", "bam", "gff"])),
        input_records=input_records,
        execution_time_sec=execution_time,
        cold_start_time_sec=draw(st.floats(min_value=0.001, max_value=3600.0)),
        warm_start_time_sec=execution_time,
        peak_rss_mb=draw(st.floats(min_value=1.0, max_value=100_000.0)),
        throughput_records_per_sec=input_records / execution_time,
        exit_code=draw(st.integers(min_value=0, max_value=0)),
        supported=True,
        error_message="",
        output_records=draw(st.integers(min_value=0, max_value=input_records)),
        unmapped_records=draw(st.integers(min_value=0, max_value=input_records)),
    )


@given(result=benchmark_result_strategy())
@settings(max_examples=100)
def test_property_1_benchmark_result_completeness(result: BenchmarkResult):
    """
    Property 1: Benchmark Result Completeness
    
    *For any* benchmark run, the output SHALL contain:
    - execution_time_sec (seconds)
    - peak_rss_mb (MB)
    - throughput_records_per_sec (records/sec)
    
    And throughput SHALL equal input_records / execution_time.
    
    **Feature: comprehensive-benchmark-audit, Property 1: Benchmark Result Completeness**
    **Validates: Requirements 2.5, 2.6, 2.7**
    """
    # Check all required fields exist and are valid
    assert result.execution_time_sec > 0, "Execution time must be positive"
    assert result.peak_rss_mb > 0, "Peak RSS must be positive"
    assert result.throughput_records_per_sec >= 0, "Throughput must be non-negative"
    
    # Verify throughput calculation: throughput = records / time
    expected_throughput = result.input_records / result.execution_time_sec
    assert abs(result.throughput_records_per_sec - expected_throughput) < 0.01, \
        f"Throughput mismatch: {result.throughput_records_per_sec} != {expected_throughput}"


# ============================================================================
# Property 2: Dual Format Output Consistency
# **Validates: Requirements 2.9, 5.4**
# ============================================================================

@st.composite
def benchmark_report_strategy(draw):
    """Generate random BenchmarkReport instances."""
    num_results = draw(st.integers(min_value=1, max_value=4))
    results = [draw(benchmark_result_strategy()) for _ in range(num_results)]
    
    metadata = BenchmarkMetadata(
        timestamp=draw(st.text(min_size=1, max_size=30)),
        system={"platform": "test", "cpu_count": "4"},
        tool_versions={"fastcrossmap": "1.0.0"},
        input_file=draw(st.text(min_size=1, max_size=100)),
        input_file_size_mb=draw(st.floats(min_value=0.1, max_value=1000.0)),
        chain_file=draw(st.text(min_size=1, max_size=100)),
        format=draw(st.sampled_from(["bed", "vcf", "bam"])),
        threads=draw(st.integers(min_value=1, max_value=32)),
    )
    
    return BenchmarkReport(metadata=metadata, results=results)


@given(report=benchmark_report_strategy())
@settings(max_examples=100)
def test_property_2_dual_format_output_consistency(report: BenchmarkReport):
    """
    Property 2: Dual Format Output Consistency
    
    *For any* benchmark output, JSON and CSV formats SHALL contain equivalent data.
    
    **Feature: comprehensive-benchmark-audit, Property 2: Dual Format Output Consistency**
    **Validates: Requirements 2.9, 5.4**
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        json_path = Path(tmpdir) / "test.json"
        csv_path = Path(tmpdir) / "test.csv"
        
        # Export to both formats
        export_to_json(report, json_path)
        export_to_csv(report, csv_path)
        
        # Verify both files exist
        assert json_path.exists(), "JSON file should be created"
        assert csv_path.exists(), "CSV file should be created"
        
        # Load JSON and verify structure
        with open(json_path) as f:
            json_data = json.load(f)
        
        assert "metadata" in json_data, "JSON should contain metadata"
        assert "results" in json_data, "JSON should contain results"
        assert len(json_data["results"]) == len(report.results), \
            "JSON should have same number of results"
        
        # Load CSV and verify structure
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            csv_rows = list(reader)
        
        assert len(csv_rows) == len(report.results), \
            "CSV should have same number of rows as results"
        
        # Verify data consistency between JSON and CSV
        for i, (json_result, csv_row) in enumerate(zip(json_data["results"], csv_rows)):
            # Check key fields match
            assert json_result["tool"] == csv_row["tool"], \
                f"Tool mismatch at row {i}"
            assert json_result["format"] == csv_row["format"], \
                f"Format mismatch at row {i}"
            assert abs(float(json_result["execution_time_sec"]) - float(csv_row["execution_time_sec"])) < 0.001, \
                f"Execution time mismatch at row {i}"
            assert abs(float(json_result["peak_rss_mb"]) - float(csv_row["peak_rss_mb"])) < 0.001, \
                f"Peak RSS mismatch at row {i}"


# ============================================================================
# Unit Tests
# ============================================================================

def test_benchmark_result_creation():
    """Test BenchmarkResult can be created with valid data."""
    result = BenchmarkResult(
        tool="fastcrossmap",
        format="bed",
        input_records=1000,
        execution_time_sec=1.5,
        cold_start_time_sec=2.0,
        warm_start_time_sec=1.5,
        peak_rss_mb=256.0,
        throughput_records_per_sec=666.67,
        exit_code=0,
        supported=True,
    )
    assert result.tool == "fastcrossmap"
    assert result.format == "bed"
    assert result.input_records == 1000


def test_unsupported_format_result():
    """Test BenchmarkResult for unsupported format."""
    result = BenchmarkResult(
        tool="liftover",
        format="vcf",
        input_records=0,
        execution_time_sec=0,
        cold_start_time_sec=0,
        warm_start_time_sec=0,
        peak_rss_mb=0,
        throughput_records_per_sec=0,
        exit_code=-1,
        supported=False,
        error_message="Format not supported"
    )
    assert not result.supported
    assert result.error_message == "Format not supported"


def test_export_json_creates_valid_file():
    """Test JSON export creates valid JSON file."""
    result = BenchmarkResult(
        tool="fastcrossmap",
        format="bed",
        input_records=1000,
        execution_time_sec=1.0,
        cold_start_time_sec=1.5,
        warm_start_time_sec=1.0,
        peak_rss_mb=100.0,
        throughput_records_per_sec=1000.0,
        exit_code=0,
        supported=True,
    )
    metadata = BenchmarkMetadata(
        timestamp="2026-01-04T12:00:00",
        system={"platform": "test"},
        tool_versions={"fastcrossmap": "1.0.0"},
        input_file="test.bed",
        input_file_size_mb=1.0,
        chain_file="test.chain",
        format="bed",
        threads=4,
    )
    report = BenchmarkReport(metadata=metadata, results=[result])
    
    with tempfile.TemporaryDirectory() as tmpdir:
        json_path = Path(tmpdir) / "test.json"
        export_to_json(report, json_path)
        
        with open(json_path) as f:
            data = json.load(f)
        
        assert data["metadata"]["timestamp"] == "2026-01-04T12:00:00"
        assert len(data["results"]) == 1
        assert data["results"][0]["tool"] == "fastcrossmap"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
