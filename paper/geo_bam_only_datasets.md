# Public Genomics Datasets Distributed as BAM (No FASTQ Available or Practical)

> Every row below was checked against a primary source (portal API, FTP README, or GEO
> SOFT record); the URL and the verbatim quote are recorded under *Verification record*.
> Rows whose claim the primary source did not support have been removed — what remains
> is the subset that holds up.
>
> The pattern across these rows is not "no FASTQ": for each resource the raw reads **do**
> exist. What varies is whether they are openly downloadable, whether a controlled-access
> application is required, and whether the BAM is the artifact the project actually
> curates. The surviving rows are **BAM-primary** or **controlled-access**, not
> "BAM-only".

## Summary table

| # | Project | Data Type | Build | Scale | Claim | Verdict |
|---|---------|-----------|-------|-------|-------|---------|
| 1 | **ENCODE** (hg19 experiments) | ChIP-seq, DNase-seq, RNA-seq, ATAC-seq | hg19 (older) / GRCh38 | 28,642 experiments; 9,018 at hg19 | "portal does host raw FASTQ for most experiments; BAM-primary, not BAM-only" | ✅ **CONFIRMED** — 164,766 FASTQ files; 26,882/28,642 (~94%) of experiments have FASTQ. BAM-primary, *not* BAM-only. |
| 2 | **GTEx** (dbGaP phs000424) | RNA-seq, WGS, WES | v7 = GRCh37/hg19; v8 = GRCh38 | ~17,000 RNA-seq samples, 948 donors | "v6/v7 hg19 BAMs via dbGaP; raw reads behind a controlled-access application" | ✅ **CONFIRMED** — v7 is GRCh37/hg19 bound to phs000424.v7.p2 (controlled). Caveat: GTEx names *BAM*, not FASTQ. |
| 3 | **UK10K** | WGS, WXS | GRCh37 (`human_g1k_v37`, `assembly=b37`) | ~10,000 individuals | "BAM files distributed via EGA under controlled access; GRCh37" | ✅ **CONFIRMED** — EGA managed access; build b37. One caveat: an EGA metadata field returned `access_type: public` for one dataset, contradicting the site's "managed access". |

## Verification record

Each block gives the exact URL fetched, a verbatim quote, what it establishes, and the
confidence. Quotes and counts marked *(independently re-checked)* were reproduced a
second time outside the verification run.

### 1. ENCODE — ✅ CONFIRMED

- **URL:** `https://www.encodeproject.org/search/?type=File&file_format=fastq&format=json&limit=0` *(independently re-checked)*
- **Quote:** `{"title": "Search", "total": 164766}` — 164,766 FASTQ files. Experiments: `total: 28642` (all), of which 26,882 have FASTQ (~94%).
- **Establishes:** the ENCODE portal does host raw FASTQ at scale, so the row is correctly labelled *BAM-primary, not BAM-only*. No correction needed.
- **Confidence:** high.
- **Note:** the WAF returns HTTP 405 "Human Verification" to browser-spoofed User-Agents; the API answers curl's default UA.

### 2. GTEx — ✅ CONFIRMED

- **URL:** `https://gtexportal.org/api/v2/metadata/dataset`; dbGaP `https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424`
- **Quote:** `{"datasetId":"gtex_v7", "dbgapId":"phs000424.v7.p2", "genomeBuild":"GRCh37/hg19"}` and `{"datasetId":"gtex_v8", "genomeBuild":"GRCh38/hg38"}`.
- **Quote:** v7 fileset "Protected Data and Raw Data" — *"Due to the nature of our donor consent agreement, raw data and attributes which might be used to identify the donors are not publicly available on the GTEx Portal. … Data available include: BAM files for RNA-Seq, Whole Exome Seq, and Whole Genome Seq …"*
- **Establishes:** v7 is hg19/GRCh37 and bound to controlled-access phs000424.v7.p2; the BAMs are in the non-public protected set. The claim holds.
- **Confidence:** high for build/accession/controlled status; medium for the exact wording "FASTQ inside dbGaP" — GTEx names *BAM*, and the raw reads are the SRA deposit behind dbGaP.
- **Could not verify:** whether dbGaP serves FASTQ directly vs. only via its linked SRA deposit (the file-type table is behind the Authorized Access login).

### 3. UK10K — ✅ CONFIRMED

- **URL:** `https://www.uk10k.org/data_access.html`; EGA `https://ega-archive.org/studies/EGAS00001000090`
- **Quote:** *"primary sequence data passing QC will be submitted to the European Genome-phenome Archive (EGA) within 60 days of generation; … all sequence/genetic variation data will be made available under a **managed access** mechanism."*
- **Quote:** Sanger VCF header — `##contig=<ID=1,assembly=b37,length=249250621>` and `##reference=…human_g1k_v37.fasta.gz`.
- **Establishes:** both raw sequences and variant calls go to EGA under managed (controlled) access; build is GRCh37, so liftover is required for GRCh38. `www.uk10k.org` is still live.
- **Confidence:** high.
- **Caveat / could not verify:** the EGA public metadata API returned `"access_type":"public"` for `EGAD00001000740` (vs `"controlled"` for `EGAD00001000195`), which contradicts the site's "managed access" wording. The API field appears unreliable, but this was not resolved from a primary source.

## Why BAM Liftover Matters

1. **Access, not absence**: for these resources the raw reads sit behind a controlled-access
   application (dbGaP / EGA), so re-aligning from FASTQ means a fresh data-access round,
   while the BAM is the curated artifact that is available once access is granted.
2. **Re-alignment is impractical**: re-aligning petabytes of data from scratch is computationally prohibitive; BAM-level coordinate liftover is orders of magnitude faster.
3. **Cross-study integration**: meta-analyses combining datasets from different genome builds (e.g. hg19 GTEx/UK10K against newer GRCh38 studies) require coordinate conversion of BAM files.
4. **Ongoing need**: as the field transitions from GRCh37/hg19 to GRCh38/hg38, large numbers of existing BAM files need coordinate conversion.
