# sequencing_comparison

## Setup

1. Create python 3 virtual environment: `virtualenv -p python3.6 venv`
2. Activate virtual environment: `source venv/bin/activate`
3. Install python packagies: `pip install -r requirements.txt`


## Step 1: prepare intervals covering CDS (protein coding exons and some padding around them)

1. cd into `intervals` directry
2. download latest GENCODE GTF file:
  
   `wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz`

3. create intervals:
  
   `python make_intervals.py -g gencode.v31.annotation.gtf.gz -o cds_intervals`
  
   Intervals will be saved into `cds_intervals` directory. By default, chromosome names have `chr` prefix. If your reference genome doesn't have this prefix, then specify `--no-chr-prefix` option when running `make_intervals.py`.
   
## Step 2: generate DP information
