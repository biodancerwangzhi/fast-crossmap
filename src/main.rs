//! FastCrossMap CLI entry point
//!
//! High-performance genome coordinate liftover tool compatible with CrossMap.

use clap::{Parser, Subcommand, ValueEnum};
use fast_crossmap::core::{ChainIndex, CoordinateMapper, ChromStyle, CompatMode};
use fast_crossmap::formats;
use std::path::PathBuf;
use std::time::Instant;

/// Compatibility mode for CrossMap behavior (CLI enum)
#[derive(Clone, Copy, Debug, Default, ValueEnum)]
pub enum CompatModeArg {
    /// Default mode: use FastCrossMap's improved logic
    #[default]
    #[value(name = "improved")]
    Improved,
    /// Strict mode: exactly match CrossMap behavior (including edge cases)
    #[value(name = "strict")]
    Strict,
}

impl From<CompatModeArg> for CompatMode {
    fn from(arg: CompatModeArg) -> Self {
        match arg {
            CompatModeArg::Improved => CompatMode::Improved,
            CompatModeArg::Strict => CompatMode::Strict,
        }
    }
}

#[derive(Parser)]
#[command(name = "fast-crossmap")]
#[command(about = "High-performance genome coordinate liftover tool")]
#[command(version)]
#[command(author = "FastCrossMap Contributors")]
struct Cli {
    /// Compatibility mode: 'strict' for CrossMap-identical output, 'improved' for enhanced logic
    #[arg(long = "compat-mode", global = true, default_value = "improved")]
    compat_mode: CompatModeArg,
    
    #[command(subcommand)]
    command: Commands,
}

#[derive(Clone, Copy, ValueEnum)]
enum ChromStyleArg {
    /// Keep chromosome names as-is
    #[value(name = "a")]
    AsIs,
    /// Use short names (1, 2, X)
    #[value(name = "s")]
    Short,
    /// Use long names (chr1, chr2, chrX)
    #[value(name = "l")]
    Long,
}

impl From<ChromStyleArg> for ChromStyle {
    fn from(arg: ChromStyleArg) -> Self {
        match arg {
            ChromStyleArg::AsIs => ChromStyle::AsIs,
            ChromStyleArg::Short => ChromStyle::Short,
            ChromStyleArg::Long => ChromStyle::Long,
        }
    }
}

#[derive(Subcommand)]
enum Commands {
    /// Convert BED format file
    Bed {
        /// Chain file for coordinate conversion
        chain: PathBuf,
        /// Input BED file
        input: PathBuf,
        /// Output file (optional, stdout if not specified)
        output: Option<PathBuf>,
        /// Number of threads (default: number of CPUs)
        #[arg(short = 't', long, default_value = "1")]
        threads: usize,
        /// Chromosome ID style: a(as-is), s(short), l(long)
        #[arg(long = "chromid", default_value = "a")]
        chrom_style: ChromStyleArg,
    },
    /// Convert VCF format file
    Vcf {
        /// Chain file for coordinate conversion
        chain: PathBuf,
        /// Input VCF file
        input: PathBuf,
        /// Output file (optional, stdout if not specified)
        output: Option<PathBuf>,
        /// Number of threads (default: number of CPUs)
        #[arg(short = 't', long, default_value = "1")]
        threads: usize,
        /// Chromosome ID style: a(as-is), s(short), l(long)
        #[arg(long = "chromid", default_value = "a")]
        chrom_style: ChromStyleArg,
    },
    /// Convert GFF/GTF format file
    Gff {
        /// Chain file for coordinate conversion
        chain: PathBuf,
        /// Input GFF/GTF file
        input: PathBuf,
        /// Output file (optional, stdout if not specified)
        output: Option<PathBuf>,
        /// Number of threads (default: number of CPUs)
        #[arg(short = 't', long, default_value = "1")]
        threads: usize,
        /// Chromosome ID style: a(as-is), s(short), l(long)
        #[arg(long = "chromid", default_value = "a")]
        chrom_style: ChromStyleArg,
    },
    /// Convert GVCF format file
    Gvcf {
        /// Chain file for coordinate conversion
        chain: PathBuf,
        /// Input GVCF file
        input: PathBuf,
        /// Output file (optional, stdout if not specified)
        output: Option<PathBuf>,
        /// Reference genome FASTA file (optional)
        #[arg(short = 'r', long)]
        refgenome: Option<PathBuf>,
        /// Compress output
        #[arg(short = 'c', long)]
        compress: bool,
        /// Number of threads (default: number of CPUs)
        #[arg(short = 't', long, default_value = "1")]
        threads: usize,
        /// Chromosome ID style: a(as-is), s(short), l(long)
        #[arg(long = "chromid", default_value = "a")]
        chrom_style: ChromStyleArg,
    },
    /// Convert MAF format file
    Maf {
        /// Chain file for coordinate conversion
        chain: PathBuf,
        /// Input MAF file
        input: PathBuf,
        /// Output file (optional, stdout if not specified)
        output: Option<PathBuf>,
        /// Reference genome FASTA file (optional)
        #[arg(short = 'r', long)]
        refgenome: Option<PathBuf>,
        /// Target genome build name (e.g., GRCh38)
        #[arg(short = 'b', long, default_value = "GRCh38")]
        build: String,
        /// Chromosome ID style: a(as-is), s(short), l(long)
        #[arg(long = "chromid", default_value = "a")]
        chrom_style: ChromStyleArg,
    },
    /// Convert Wiggle/bedGraph format file
    Wig {
        /// Chain file for coordinate conversion
        chain: PathBuf,
        /// Input Wiggle/bedGraph file
        input: PathBuf,
        /// Output file (optional, stdout if not specified)
        output: Option<PathBuf>,
        /// Chromosome ID style: a(as-is), s(short), l(long)
        #[arg(long = "chromid", default_value = "a")]
        chrom_style: ChromStyleArg,
    },
    /// Convert BAM/SAM/CRAM format file
    #[cfg(feature = "bam")]
    Bam {
        /// Chain file for coordinate conversion
        chain: PathBuf,
        /// Input BAM/SAM/CRAM file
        input: PathBuf,
        /// Output BAM file
        output: PathBuf,
        /// Number of threads for parallel I/O
        #[arg(short = 't', long, default_value = "1")]
        threads: usize,
        /// Chromosome ID style: a(as-is), s(short), l(long)
        #[arg(long = "chromid", default_value = "a")]
        chrom_style: ChromStyleArg,
    },
    /// Convert large genomic regions (partial mapping allowed)
    Region {
        /// Chain file for coordinate conversion
        chain: PathBuf,
        /// Input BED file with regions
        input: PathBuf,
        /// Output file (optional, stdout if not specified)
        output: Option<PathBuf>,
        /// Minimum mapping ratio (default: 0.85)
        #[arg(short = 'r', long, default_value = "0.85")]
        ratio: f64,
        /// Chromosome ID style: a(as-is), s(short), l(long)
        #[arg(long = "chromid", default_value = "a")]
        chrom_style: ChromStyleArg,
    },
    /// Convert BigWig format file
    Bigwig {
        /// Chain file for coordinate conversion
        chain: PathBuf,
        /// Input BigWig file
        input: PathBuf,
        /// Output file prefix (will create .bgr file)
        output: Option<PathBuf>,
        /// Chromosome ID style: a(as-is), s(short), l(long)
        #[arg(long = "chromid", default_value = "a")]
        chrom_style: ChromStyleArg,
    },
}


fn load_chain(chain_path: &PathBuf, chrom_style: ChromStyleArg, compat_mode: CompatModeArg) -> anyhow::Result<CoordinateMapper> {
    let start = Instant::now();
    eprintln!("Loading chain file: {:?}", chain_path);
    
    let index = ChainIndex::from_chain_file(chain_path)
        .map_err(|e| anyhow::anyhow!("Failed to load chain file: {}", e))?;
    
    let mapper = CoordinateMapper::with_compat_mode(index, chrom_style.into(), compat_mode.into());
    eprintln!("Chain file loaded in {:.2}s", start.elapsed().as_secs_f64());
    
    Ok(mapper)
}

fn main() -> anyhow::Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    let start = Instant::now();
    
    // Log compatibility mode
    match cli.compat_mode {
        CompatModeArg::Strict => eprintln!("Compatibility mode: strict (CrossMap-identical)"),
        CompatModeArg::Improved => {} // Don't log for default mode
    }

    match cli.command {
        Commands::Bed { chain, input, output, threads, chrom_style } => {
            let mapper = load_chain(&chain, chrom_style, cli.compat_mode)?;
            let output_path = output.unwrap_or_else(|| PathBuf::from("output.bed"));
            let unmap_path = output_path.with_extension("bed.unmap");
            
            eprintln!("Converting BED file: {:?} -> {:?}", input, output_path);
            let stats = formats::convert_bed(&input, &output_path, &unmap_path, &mapper, threads)?;
            
            eprintln!("\n=== Conversion Statistics ===");
            eprintln!("Total records:   {}", stats.total);
            eprintln!("Successful:      {}", stats.success);
            eprintln!("Failed:          {}", stats.failed);
            eprintln!("Time elapsed:    {:.2}s", start.elapsed().as_secs_f64());
        }
        
        Commands::Vcf { chain, input, output, threads, chrom_style } => {
            let mapper = load_chain(&chain, chrom_style, cli.compat_mode)?;
            let output_path = output.unwrap_or_else(|| PathBuf::from("output.vcf"));
            
            eprintln!("Converting VCF file: {:?} -> {:?}", input, output_path);
            let stats = formats::convert_vcf(&input, &output_path, &mapper, None::<&PathBuf>, false, threads)?;
            
            eprintln!("\n=== Conversion Statistics ===");
            eprintln!("Total records:   {}", stats.total);
            eprintln!("Successful:      {}", stats.success);
            eprintln!("Failed:          {}", stats.failed);
            eprintln!("Time elapsed:    {:.2}s", start.elapsed().as_secs_f64());
        }
        
        Commands::Gff { chain, input, output, threads, chrom_style } => {
            let mapper = load_chain(&chain, chrom_style, cli.compat_mode)?;
            let output_path = output.unwrap_or_else(|| PathBuf::from("output.gff"));
            
            eprintln!("Converting GFF file: {:?} -> {:?}", input, output_path);
            let stats = formats::convert_gff(&input, &output_path, &mapper, threads)?;
            
            eprintln!("\n=== Conversion Statistics ===");
            eprintln!("Total records:   {}", stats.total);
            eprintln!("Successful:      {}", stats.success);
            eprintln!("Failed:          {}", stats.failed);
            eprintln!("Time elapsed:    {:.2}s", start.elapsed().as_secs_f64());
        }
        
        Commands::Gvcf { chain, input, output, refgenome, compress, threads, chrom_style } => {
            let mapper = load_chain(&chain, chrom_style, cli.compat_mode)?;
            let output_path = output.unwrap_or_else(|| PathBuf::from("output.gvcf"));
            
            eprintln!("Converting GVCF file: {:?} -> {:?}", input, output_path);
            let stats = formats::convert_gvcf(
                &input, &output_path, &mapper, 
                refgenome.as_ref(), compress, threads
            )?;
            
            eprintln!("\n=== Conversion Statistics ===");
            eprintln!("Total records:   {}", stats.total);
            eprintln!("Successful:      {}", stats.success);
            eprintln!("Failed:          {}", stats.failed);
            eprintln!("Time elapsed:    {:.2}s", start.elapsed().as_secs_f64());
        }
        
        Commands::Maf { chain, input, output, refgenome, build, chrom_style } => {
            let mapper = load_chain(&chain, chrom_style, cli.compat_mode)?;
            let output_path = output.unwrap_or_else(|| PathBuf::from("output.maf"));
            
            eprintln!("Converting MAF file: {:?} -> {:?}", input, output_path);
            let stats = formats::convert_maf(
                &input, &output_path, &mapper, 
                refgenome.as_ref(), &build
            )?;
            
            eprintln!("\n=== Conversion Statistics ===");
            eprintln!("Total records:   {}", stats.total);
            eprintln!("Successful:      {}", stats.success);
            eprintln!("Failed:          {}", stats.failed);
            eprintln!("Time elapsed:    {:.2}s", start.elapsed().as_secs_f64());
        }
        
        Commands::Wig { chain, input, output, chrom_style } => {
            let mapper = load_chain(&chain, chrom_style, cli.compat_mode)?;
            let output_path = output.unwrap_or_else(|| PathBuf::from("output.bedGraph"));
            
            eprintln!("Converting Wiggle file: {:?} -> {:?}", input, output_path);
            let stats = formats::convert_wig(&input, &output_path, &mapper)?;
            
            eprintln!("\n=== Conversion Statistics ===");
            eprintln!("Total records:   {}", stats.total);
            eprintln!("Successful:      {}", stats.success);
            eprintln!("Failed:          {}", stats.failed);
            eprintln!("Merged:          {}", stats.merged);
            eprintln!("Time elapsed:    {:.2}s", start.elapsed().as_secs_f64());
        }
        
        #[cfg(feature = "bam")]
        Commands::Bam { chain, input, output, threads, chrom_style } => {
            let mapper = load_chain(&chain, chrom_style, cli.compat_mode)?;
            
            eprintln!("Converting BAM file: {:?} -> {:?}", input, output);
            let stats = formats::convert_bam(&input, &output, &mapper, threads)?;
            
            eprintln!("\n=== Conversion Statistics ===");
            eprintln!("Total records:   {}", stats.total);
            eprintln!("Mapped:          {}", stats.mapped);
            eprintln!("Unmapped:        {}", stats.unmapped);
            eprintln!("Failed:          {}", stats.failed);
            eprintln!("Paired:          {}", stats.paired);
            eprintln!("Single:          {}", stats.single);
            eprintln!("Time elapsed:    {:.2}s", start.elapsed().as_secs_f64());
        }
        
        Commands::Region { chain, input, output, ratio, chrom_style } => {
            let mapper = load_chain(&chain, chrom_style, cli.compat_mode)?;
            let output_path = output.unwrap_or_else(|| PathBuf::from("output.bed"));
            
            eprintln!("Converting Region file: {:?} -> {:?} (min_ratio={})", input, output_path, ratio);
            let stats = formats::convert_region(&input, &output_path, &mapper, ratio)?;
            
            eprintln!("\n=== Conversion Statistics ===");
            eprintln!("Total records:   {}", stats.total);
            eprintln!("Successful:      {}", stats.success);
            eprintln!("Failed:          {}", stats.failed);
            eprintln!("  - Unmapped:    {}", stats.unmapped);
            eprintln!("  - CrossChrom:  {}", stats.cross_chrom);
            eprintln!("  - LowRatio:    {}", stats.low_ratio);
            eprintln!("Time elapsed:    {:.2}s", start.elapsed().as_secs_f64());
        }
        
        Commands::Bigwig { chain, input, output, chrom_style } => {
            let mapper = load_chain(&chain, chrom_style, cli.compat_mode)?;
            let output_path = output.unwrap_or_else(|| PathBuf::from("output"));
            
            eprintln!("Converting BigWig file: {:?} -> {:?}", input, output_path);
            let stats = formats::convert_bigwig(&input, &output_path, &mapper, false)?;
            
            eprintln!("\n=== Conversion Statistics ===");
            eprintln!("Total records:   {}", stats.total);
            eprintln!("Successful:      {}", stats.success);
            eprintln!("Failed:          {}", stats.failed);
            eprintln!("Merged:          {}", stats.merged);
            eprintln!("Time elapsed:    {:.2}s", start.elapsed().as_secs_f64());
        }
    }

    Ok(())
}
