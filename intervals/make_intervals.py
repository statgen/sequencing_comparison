import os
import argparse
import gzip
from intervaltree import IntervalTree


argparser = argparse.ArgumentParser(description = 'Creates file with coding (CDS) intervals based on GENCODE GTF.')
argparser.add_argument('-g', '--gencode', metavar = 'file', dest = 'in_gencode_file', required = True, help = 'Input GENCODE file in GTF format.')
argparser.add_argument('-o', '--output', metavar = 'directory', dest = 'out_directory', required = True, help = 'Output directory name.')
argparser.add_argument('-p', '--padding', metavar = 'number', dest = 'in_padding', required = False, default = 1000, type = int, help = 'Padding around genes in base-pairs. Default: 1000.')
argparser.add_argument('--no-chr-prefix', dest = 'no_chr_prefix', action = 'store_true', default = False, help = 'Remove "chr" prefix.')

chromosomes = list(map(str, range(1, 23))) + ['X']

def write_bucket(directory, bucket, no_chr_prefix):
    chrom = bucket['chrom']
    if no_chr_prefix and chrom.startswith('chr'):
        chrom = chrom[3:]
    elif not no_chr_prefix and not chrom.startswith('chr'):
        chrom = 'chr' + chrom
    filename = os.path.join(directory, f'{chrom}-{bucket["start"]}-{bucket["stop"]}.list')
    with open(filename, 'w') as ofile:
        for bucket_interval in bucket['intervals']:
            ofile.write(f'{chrom}\t{bucket_interval.begin}\t{bucket_interval.end - 1}\n') # end - 1 because end does not include upper bound

if __name__ == '__main__':
    args = argparser.parse_args()

    intervals_by_chrom = dict()
    with gzip.open(args.in_gencode_file, 'rt') as ifile:
        for line in ifile:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            feature_type = fields[2]
            # start codon is already included to CDS
            # stop codon is not included to CDS by GENCODE, but lets not consider it at this point, because it will be included with padding
            if feature_type != 'CDS':
                continue
            chrom = fields[0]
            start_bp, stop_bp = map(int, fields[3:5])
            attributes = dict(map(lambda x: x.strip('"'), x.strip().split()) for x in fields[8].split(';') if x.strip() != '')
            if attributes['gene_type'] != 'protein_coding' or attributes['transcript_type'] != 'protein_coding':
                continue
            # skip if the coding region end or start could not be confirmed
            if any(x in ['cds_end_NF', 'cds_start_NF', 'mRNA_end_NF', 'mRNA_start_NF'] for x in attributes.get('tag', [])):
                continue
            intervals = intervals_by_chrom.setdefault(chrom, IntervalTree())
            intervals.addi(max(0, start_bp - args.in_padding), stop_bp + args.in_padding + 1) # does not include upper bound

    os.mkdir(args.out_directory)

    for chrom in chromosomes:
        intervals = intervals_by_chrom.get(chrom, None)
        if intervals is None:
            chrom = 'chr' + chrom
            intervals = intervals_by_chrom.get(chrom, None)
            if intervals is None:
                continue
        intervals.merge_overlaps()
        bucket = {}
        for interval in sorted(intervals):
            if bucket:
                gap_bp = interval.begin - bucket['stop']
                if gap_bp <= 0:
                    raise Exception("Overlapping intervals detected.")
                elif gap_bp < 100000 and (interval.end - bucket['start']) < 500000: # if gap is small and bucket is not larger than 500kb then add interval to the bucket
                    bucket['stop'] = interval.end - 1 # end doesn't include upper bound
                    bucket['intervals'].append(interval)
                    continue
                else:
                    write_bucket(args.out_directory, bucket, args.no_chr_prefix)
                    bucket = {}
            bucket['chrom'] = chrom
            bucket['start'] = interval.begin
            bucket['stop'] = interval.end - 1 # end doesn't include upper bound
            bucket['intervals'] = [ interval ]
        if bucket:
            write_bucket(args.out_directory, bucket, args.no_chr_prefix)
