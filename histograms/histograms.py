#!/usr/bin/env python3


import argparse
import gzip
import pysam
from intervaltree import IntervalTree
from collections import Counter
import json


argparser = argparse.ArgumentParser(description = 'Histograms of sequencing depth across CDS.')
argparser.add_argument('-g', '--gencode', metavar = 'file', dest = 'in_gencode_file', required = True, help = 'Input GENCODE file in GTF format.')
argparser.add_argument('-d', '--depth', metavar = 'file', dest = 'in_vcf', required = True, help = 'Input VCF/BCF with DP FORMAT field.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_file', required = True, help = 'Summary file in compressed JSON format.')


if __name__ == '__main__':
    args = argparser.parse_args()

    genes_by_chrom = dict()
    with gzip.open(args.in_gencode_file, 'rt') as ifile:
        for line in ifile:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            feature_type = fields[2]
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
            genes = genes_by_chrom.setdefault(chrom, dict())
            gene_info = genes.setdefault(attributes['gene_id'], { 'id': attributes['gene_id'], 'name': attributes['gene_name'], 'cds': IntervalTree() })
            gene_info['cds'].addi(start_bp, stop_bp + 1)

    with pysam.VariantFile(args.in_vcf, 'r') as ifile, gzip.open(args.out_file, 'wt') as ofile:
        contigs = list(ifile.header.contigs)
        jsons = []
        for chrom, genes in genes_by_chrom.items():
            if chrom not in contigs:
                if chrom.startswith('chr'):
                    if chrom[3:] not in contigs:
                        continue
                elif 'chr' + chrom not in contigs:
                    continue
            cds_all = IntervalTree()
            cds_gene = IntervalTree()
            histogram_all = Counter()
            histogram_gene = dict()
            for gene_name, gene_info in genes.items():
                gene_info['cds'].merge_overlaps()
                histogram_gene[gene_info['id']] = Counter()
                for interval in gene_info['cds']:
                    cds_all.addi(interval.begin, interval.end)
                    cds_gene.addi(interval.begin, interval.end, gene_info['id'])
            cds_all.merge_overlaps()
            for interval in cds_all:
                records = list(ifile.fetch(chrom, interval.begin - 1, interval.end - 1))
                if not records: # if nothing fetched, try adding or removing 'chr' prefix
                    if chrom.startswith('chr'):
                        records = ifile.fetch(chrom[3:], interval.begin - 1, interval.end - 1)
                    else:
                        records = ifile.fetch('chr' + chrom, interval.begin - 1, interval.end - 1)
                for record in records:
                    overlapping_genes = [x.data for x in cds_gene.at(record.pos)]
                    dps = [fmt['DP'] for _, fmt in record.samples.items()]
                    histogram_all.update(dps)
                    for overlapping_gene in overlapping_genes:
                        histogram_gene[overlapping_gene].update(dps)
            if len(histogram_all) > 0:
                jsons.append({'chrom': chrom, 'type': 'all', 'histogram': histogram_all})
                for gene_id, gene_histogram in histogram_gene.items():
                    jsons.append({'chrom': chrom, 'type': 'gene', 'gene_id': gene_id, 'histogram': gene_histogram})
        json.dump(jsons, ofile)
