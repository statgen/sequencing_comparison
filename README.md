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

1. cd into `coverage` directory
2. edit `nextflow.config` file:
   
   `bams_list_path` -- Same file as in coverage computation: points to the file were each line is a whitespace delimited tuple: sample name, absolute path to  the corresponding BAM/CRAM file.
   
   `vcfs` -- path to VCF or BCF files with genotype information. Files must include both PASS and QC failed variants.
   
   `reference_path` -- points to FASTA file (`*.fa`) with genome reference. File must be indexed (ie. the corresponding `*.fai` should in the same directory)
   
   `bcftools` -- path to `bcftools` executable
   
   `tabix` -- path to `tabix` executable
   
   Edit other options related to SLURM or local execution as needed.
   
3. Run `nextlow run Subset.nf`. Preferrably run from `tmux` session. When crashed (e.g. SLURM node failure) use `nextflow run Subset.nf -resume`.

4. The final BCF files with normalized variants are located in `results/` folder.
