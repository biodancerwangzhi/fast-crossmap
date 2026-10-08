# Public Genomics Datasets Distributed as BAM (No FASTQ Available or Practical)

This table lists major public genomics resources where data is primarily or exclusively
distributed as BAM files, making BAM-level coordinate liftover tools essential.

## Large-Scale Projects with BAM-Only or BAM-Primary Distribution

| # | Project / Accession | Data Type | Genome Build | Scale | BAM-Only Reason | Reference |
|---|---------------------|-----------|-------------|-------|----------------|-----------|
| 1 | **TCGA (GDC Legacy Archive)** | WGS, WXS, RNA-seq | hg19 (legacy) / hg38 (harmonized) | ~11,000 cases, 33 tumor types, ~2.5 PB total | Legacy BAMs aligned to hg19; FASTQ not provided for legacy data. GDC re-aligned to hg38 but many users still work with legacy hg19 BAMs. | [GDC Legacy vs Harmonized](https://gdc.cancer.gov/about-data/publications/HG38QC); [Cell Systems 2019](https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30201-7) |
| 2 | **GTEx (dbGaP phs000424)** | RNA-seq, WGS | hg19 (v6/v7), hg38 (v8+) | ~17,000 RNA-seq samples, 948 donors, 54 tissues | Earlier releases (v6, v7) distributed as hg19-aligned BAMs via dbGaP; FASTQ requires SRA conversion. Many labs still use v7 hg19 BAMs. | [GTEx Portal](https://gtexportal.org); [dbGaP phs000424](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424) |
| 3 | **1000 Genomes Phase 3** | WGS (low-coverage), Exome | hg19 (GRCh37) | 2,504 individuals, ~100 TB BAMs | Data distributed as BAMs on S3/EBI FTP; no FASTQ provided. BAMs aligned to GRCh37, need liftover for GRCh38 analyses. | [1000 Genomes S3](https://1000genomes.s3.amazonaws.com); [IGSR](https://www.internationalgenome.org/data/) |
| 4 | **ENCODE (hg19 experiments)** | ChIP-seq, DNase-seq, RNA-seq, ATAC-seq | hg19 (older) / GRCh38 (newer) | >15,000 experiments, thousands of BAM files | Many early ENCODE experiments (2012–2015) have hg19-aligned BAMs as primary processed files; some lack raw FASTQ on GEO. | [ENCODE Portal](https://www.encodeproject.org); [ENCODE guidelines, Genome Res 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3431496/) |
| 5 | **Roadmap Epigenomics** | ChIP-seq, DNase-seq, RNA-seq, WGBS | hg19 | 111 reference epigenomes, >2,800 experiments | All processed BAMs aligned to hg19; FASTQ scattered across SRA. Most users download pre-aligned BAMs from the Roadmap data portal. | [Roadmap Epigenomics](https://egg2.wustl.edu/roadmap/web_portal/) |
| 6 | **ICGC (International Cancer Genome Consortium)** | WGS, WXS | hg19 (GRCh37) | >25,000 tumor-normal pairs, 89 cancer projects | Many ICGC projects distribute aligned BAMs only; raw data restricted. BAMs aligned to GRCh37 need liftover for pan-cancer GRCh38 studies. | [ICGC Data Portal](https://dcc.icgc.org) |
| 7 | **BLUEPRINT Epigenome** | ChIP-seq, DNase-seq, RNA-seq, WGBS | hg19 (GRCh37) | >400 samples, haematopoietic cells | Processed BAMs aligned to GRCh37 are the primary download format; FASTQ not always available. | [BLUEPRINT](http://www.blueprint-epigenome.eu) |
| 8 | **TARGET (Therapeutically Applicable Research to Generate Effective Treatments)** | WGS, WXS, RNA-seq | hg19 (legacy) / hg38 (harmonized) | ~6,200 cases, pediatric cancers | Similar to TCGA: legacy BAMs in hg19 on GDC. FASTQ not provided for legacy data. | [GDC TARGET](https://portal.gdc.cancer.gov/projects?filters=%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22projects.program.name%22%2C%22value%22%3A%22TARGET%22%7D%7D) |
| 9 | **CGCI (Cancer Genome Characterization Initiative)** | WGS, WXS, RNA-seq | hg19 / hg38 | Multiple cancer types | Legacy BAMs on GDC in hg19; harmonized BAMs in hg38. No FASTQ for legacy data. | [GDC Portal](https://portal.gdc.cancer.gov) |
| 10 | **UK10K** | WGS, WXS | hg19 (GRCh37) | ~10,000 individuals | BAM files distributed via EGA; FASTQ not always included. Aligned to GRCh37, requires liftover for GRCh38. | [UK10K](https://www.uk10k.org) |

## GEO/SRA Examples with BAM-Only Submissions

| # | Accession | Data Type | Species | Samples | Note |
|---|-----------|-----------|---------|---------|------|
| 11 | **GSE118165** | scRNA-seq (10x) | Human | ~40 | BAM submitted to SRA; FASTQ reconstruction needed via bamtofastq |
| 12 | **GSE122960** | scRNA-seq (10x) | Human | 8 lungs | [Only BAM in SRA, no direct FASTQ FTP links](https://www.biostars.org/p/470724/) |
| 13 | **GSE128639** | scRNA-seq (10x) | Human | Multiple | BAM-only SRA submission; common for 10x Chromium data |
| 14 | **GSE84133** | scRNA-seq (inDrop) | Human/Mouse | Pancreatic islets | BAM files only; noted on [Biostars](https://www.biostars.org/p/9556827/) as common pattern |
| 15 | **GSE81547** | scRNA-seq | Human | Pancreas | BAM submitted, no FASTQ in SRA |

## Summary Statistics

- **Total BAM-only data volume**: >5 petabytes across these projects
- **Number of samples with hg19 BAMs needing potential liftover**: >50,000
- **Key observation**: SRA [accepts BAM as a valid submission format](https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/), and many submitters — especially for 10x Genomics single-cell data — submit only BAM files. NCBI [stores ETL (Extract-Transform-Load) data](https://helixweb.nih.gov/apps/sratoolkit.html), and original BAMs are only available from AWS/GCP cloud mirrors.

## Why BAM Liftover Matters

1. **No FASTQ fallback**: For many legacy datasets (TCGA, 1000G Phase 3, Roadmap), raw FASTQ files are not available — users must work directly with BAMs.
2. **Re-alignment is impractical**: Re-aligning petabytes of data from scratch is computationally prohibitive. BAM-level coordinate liftover is orders of magnitude faster.
3. **Cross-study integration**: Meta-analyses combining datasets from different genome builds (e.g., TCGA hg19 + newer GRCh38 studies) require coordinate conversion of BAM files.
4. **Ongoing need**: As the field transitions from GRCh37/hg19 to GRCh38/hg38, millions of existing BAM files need coordinate conversion.
