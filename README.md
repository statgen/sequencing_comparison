# sequencing_comparison

## A. Compute coverage

### Setup

1. Make sure that the latest versions of `samtools`, `bcftools`, and `tabix` are installed.
2. Download the latest Nextflow from `https://www.nextflow.io`.
3. Create python 3 virtual environment: `virtualenv -p python3.6 venv`.
4. Activate virtual environment: `source venv/bin/activate`.
5. Install python packagies: `pip install -r requirements.txt`.


### Run Step 1: prepare intervals covering CDS regions (i.e. protein coding exons and some padding around them)

1. cd into `intervals` directry
2. download latest GENCODE GTF file:
  
   `wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz`

3. create intervals:
  
   `python make_intervals.py -g gencode.v31.annotation.gtf.gz -o cds_intervals`
  
   Intervals will be saved into `cds_intervals` directory. By default, chromosome names have `chr` prefix. If your reference genome doesn't have this prefix, then specify `--no-chr-prefix` option when running `make_intervals.py`.
   
### Run Step 2: generate DP information

1. cd into `coverage` directory
2. edit `nextflow.config` file:
   
   `bams_list_path` -- points to the file were each line is a whitespace delimited tuple: sample name, absolute path to  the corresponding BAM/CRAM file.
   
   `cram` -- true if working with CRAMs
   
   `intervals` -- points to `*.list` files with intervals generated in Step 1.
   
   `reference_path` -- points to FASTA file (`*.fa`) with genome reference. File must be indexed (ie. the corresponding `*.fai` should in the same directory)
   
   `samtools` -- path to `samtools` executable
   
   `max_depth` -- maximal depth per sample (see documentation for `samtools mpileup` for details)
   
   `mpileup_subset` -- points to `mpileup_subset.py` script (prefferable absolute path)
   
   `bcftools` -- path to `bcftools` executable
   
   `tabix` -- path to `tabix` executable
   
   Edit other options related to SLURM or local execution as needed.
   
3. Run `nextlow run Coverage.nf`. Preferrably run from `tmux` session. When crashed (e.g. SLURM node failure) use `nextflow run Coverage.nf -resume`.

4. The final BCF files with DP information for each CDS base-pair and each sample are located in `results/merged` folder.

## B. Subset VCF and normalize variants

### Setup

1. Make sure that the latest versions of `bcftools` and `tabix` are installed.
2. Download the latest Nextflow from `https://www.nextflow.io`.

### Run

1. cd into `subset_normalize` directory
2. edit `nextflow.config` file:
   
   `bams_list_path` -- Same file as in coverage computation: points to the file were each line is a whitespace delimited tuple: sample name, absolute path to  the corresponding BAM/CRAM file.
   
   `vcfs` -- path to VCF or BCF files with genotype information. Files must include both PASS and QC failed variants.
   
   `reference_path` -- points to FASTA file (`*.fa`) with genome reference. File must be indexed (ie. the corresponding `*.fai` should in the same directory)
   
   `bcftools` -- path to `bcftools` executable
   
   `tabix` -- path to `tabix` executable
   
   Edit other options related to SLURM or local execution as needed.
   
3. Run `nextlow run Subset.nf`. Preferrably run from `tmux` session. When crashed (e.g. SLURM node failure) use `nextflow run Subset.nf -resume`.

4. The final BCF files with normalized variants are located in `results/` folder.

## C. Annotate variants

### Setup

1. Make sure that the latest version of `bcftools` and `tabix` are installed
2. Make sure that the `VEP` (Variant Effect Predictor) v96 or higher is installed.

### Run

1. cd into `annotate` directory
2. edit `nextflow.config` file:

   `vcfs` -- path to the VCF/BCF files generate in the previous step
   
   `vep` -- path to `VEP` executable
   
   `vep_flags` -- set any additional VEP flags if needed
   
   `bcftools` -- path to `VEP` executable
   
   `tabix` -- path to `tabix` executable

   Edit other options related to SLURM or local execution as needed.
   
3. Run `nextlow run Annotate.nf`. Preferrably run from `tmux` session. When crashed (e.g. SLURM node failure) use `nextflow run Annotate.nf -resume`.

4. The final VCF files with annotated variants are located in `results/vep` folder.

## D. Pairwise sample comparison (concordance)

1. cd into `compare` directory
2. edit `nextlow.config` file:

   `pairs_list_path` -- Mapping file with sample IDs (no header). Each line has two whitespace (or tab) delimited columns. First column stores sample ID in study 1, second column stores corresponding sample ID in study 2. See `example_pairs.list`.
   
   `study1_files_list_path` -- File with absolute paths to coverage, genotype, and annotation files generated in previous steps for study 1. The file has three whitespace (or tab) delimited columns (no header). The first column stores absolute path to coverage files from step A, the second column stores absolute path to the genotype files from step B, the third column stores absolute path to the annotation files from step C. Important: each row must store corresponding files for the same chromosome. See `study1_files.list`.
   
   `study2_files_list_path` -- File with the same structure as `study1_files_list_path`, but with absolute paths to coverage, genotype, and annotation files for study 2.
   
   `compare` -- absolute path to the`compare.py` script
   
   Edit other options related to SLURM or local execution as needed.
   
3. Run `nextlow run Compare.nf`. Preferrably run from `tmux` session. When crashed (e.g. SLURM node failure) use `nextflow run Compare.nf -resume`.

4. The final gzip compressed summary files are located in `results/` folder.

## E. Depth histograms in CDS

1. cd into `histograms` directory
2. edit `nextflow.config`
  
   `gencode_gtf_path` -- path to GENCODE GTF file. Change only if other than `v31` version needed.
  
   `coverage_files_path` -- path to VCF/BCF files generated in step A (i.e. in `sequencing_comparison/coverage/results/merged` directory).
  
   `coverage_files_index_suffix` -- change to `tbi` if VCF/BCF files in step A were indexed using TBI index (i.e. default tabix).
  
   `histograms` -- absolute path to the`histograms.py` script.
  
   Edit other options related to SLURM or local execution as needed.
  
3. Run `nextlow run Histograms.nf`. Preferrably run from `tmux` session. When crashed (e.g. SLURM node failure) use `nextflow run Histograms.nf -resume`.

4. The final gzip compressed histogram files are located in `results/` folder.
