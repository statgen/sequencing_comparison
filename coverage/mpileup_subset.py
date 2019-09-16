#!/usr/bin/env python3

import argparse
import sys
import pysam


argparser = argparse.ArgumentParser(description = 'Subsets DEPTH columns from multi-sample samtools mpileup output.')
argparser.add_argument('-p', '--pileup', metavar = 'file', dest = 'in_mpileup_file', type = argparse.FileType('rt'), required = True, help = 'Input samtools mpileup (uncompressed). Specify \'-\' when reading from stdin.')
argparser.add_argument('-s', '--samples', metavar = 'name', dest = 'in_samples', type = str, nargs = '+', required = True, help = 'List of sample names. Order must correspond to samtools mpileup input BAMs/CRAMs.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'output_file', required = True, help = 'Output file compressed using bgzip.')


if __name__ == '__main__':
    args = argparser.parse_args()
    header = pysam.VariantHeader()
    header.add_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of bases">')
    for chrom in list(map(str, range(1, 23))) + ['X', 'Y']: # add all possible chromosome names. needed by htslibs to track chromosomes in records.
        header.contigs.add(chrom)
        header.contigs.add('chr' + chrom)
    for sample in args.in_samples:
        header.add_sample(sample)
    with pysam.VariantFile(args.output_file, 'wb', header = header) as ofile:
        for line in args.in_mpileup_file:
            fields = line.rstrip().split('\t')
            record = header.new_record(contig = fields[0], start = int(fields[1]) - 1, stop = int(fields[1]), alleles = ( fields[2], '.' ), samples =  [ { 'DP': int(fields[i]) } for i in range(3, len(fields), 3) ] )
            record.info.clear()
            ofile.write(record)
