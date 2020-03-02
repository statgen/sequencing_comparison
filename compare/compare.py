#!/usr/bin/env python3

import argparse
import pysam
import statistics
import json
import gzip
import re
from intervaltree import IntervalTree
import time
from collections import OrderedDict


argparser = argparse.ArgumentParser(description = 'Compares called genotypes in two samples from different sequencing experiments')
argparser.add_argument('-s1', '--sample-study1', metavar = 'name', dest = 'sample1', required = True, help = 'Sample name from study 1.')
argparser.add_argument('-g1', '--genotypes-study1', metavar = 'file', dest = 'in_vcf1', required = True, help = 'Genotypes from study 1 in VCF/BCF.')
argparser.add_argument('-d1', '--depth-study1', metavar = 'file', dest = 'in_dp_vcf1', required = True, help = 'Depth (DP) from study 1 in VCF/BCF.')
argparser.add_argument('-a1', '--annotations-study1', metavar = 'file', dest = 'in_vep_vcf1', required = True, help = 'VEP annotations from study 1 in VCF/BCF.')
argparser.add_argument('-s2', '--sample-study2', metavar = 'name', dest = 'sample2', required = True, help = 'Sample name from study 2.')
argparser.add_argument('-g2', '--genotypes-study2', metavar = 'file', dest = 'in_vcf2', required = True, help = 'Genotypes from study 2 in VCF/BCF.')
argparser.add_argument('-d2', '--depth-study2', metavar = 'file', dest = 'in_dp_vcf2', required = True, help = 'Depth (DP) from study 2 in VCF/BCF.')
argparser.add_argument('-a2', '--annotations-study2', metavar = 'file', dest = 'in_vep_vcf2', required = True, help = 'VEP annotations from study 2 in VCF/BCF.')
argparser.add_argument('-t', '--targets', metavar = 'file', dest = 'in_targets_bed', required = True, help = 'BED file with target regions (e.g. target regions in WES).')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_file', required = True, help = 'Output file prefix.')


class Variant(object):
    __slots__ = ['is_pass', 'an', 'ac', 'gt', 'dp', 'cat']

    def __init__(self, gt, is_pass):
        self.gt = gt
        self.is_pass = is_pass
        self.an = None
        self.ac = None
        self.dp = None
        self.cat = None

    def __str__(self):
        return f'{self.is_pass}\t{self.an}\t{self.ac}\t{self.gt}\t{self.dp}\t{self.cat}'


cds_variant_types = [
        'start_lost',
        'start_retained_variant',
        'stop_retained_variant',
        'synonymous_variant',
        'missense_variant',
        'inframe_insertion',
        'inframe_deletion',
        'stop_gained',
        'stop_lost',
        'frameshift_variant',
        'splice_donor_variant',
        'splice_acceptor_variant'
        ]


alt_regex = re.compile('^[AGTC]+$')


def read_targets(in_bed, targets):
    with open(in_bed, 'rt') as ifile:
        for line in ifile:
            fields = line.rstrip().split()
            if len(fields) < 3:
                continue
            chrom = fields[0][3:] if fields[0].startswith('chr') else fields[0]
            start = int(fields[1]) + 1 # bed encodes first chromosomal position as 0
            stop = int(fields[2]) # bed stores open intervals, so no need to add 1
            chrom_targets = targets.setdefault(chrom, IntervalTree())
            chrom_targets.addi(start, stop + 1) # IntervalTree stores open end intevals, so we need to add 1 to stop.


def read_genotypes(in_vcf, sample, variants):
    with pysam.VariantFile(in_vcf) as ifile:
        vcf_samples = set((ifile.header.samples))
        if not sample in vcf_samples:
            raise Exception(f'Sample {sample} was not found in {in_vcf} file.')
        ifile.subset_samples([sample])
        for record in ifile:
            if len(record.alts) > 1: # multi-allelic variants must be split into multiple bi-allelic VCF entries
                raise Exception('Multi-allelic VCF records are not supported. Multi-allelic variants must be split into multiple bi-allelic VCF entries.')
            alt = record.alts[0]
            if not alt_regex.match(alt): # skip non-trivial variants e.g. ALT=* for spanning deletions
                continue
            gt = record.samples[sample]['GT']
            if None in gt or sum(gt) == 0:
                continue
            variant_name = (record.chrom[3:] if record.chrom.startswith('chr') else record.chrom, record.pos, record.ref, alt)
            variants[variant_name] = Variant(gt = sum(gt), is_pass = 'PASS' in list(record.filter))


def read_vep(in_vep_vcf, variants):
    with pysam.VariantFile(in_vep_vcf) as ifile:
        csq_meta = ifile.header.info.get('CSQ', None)
        if csq_meta is None:
            raise Exception('No meta-information entry about CSQ INFO field found!')
        csq_header = csq_meta.description.split(':', 1)[1].strip().split('|')
        for record in ifile.fetch():
            if len(record.alts) > 1: # multi-allelic variants must be split into multiple bi-allelic VCF entries
                raise Exception('Multi-allelic VCF records are not supported. Multi-allelic variants must be split into multiple bi-allelic VCF entries.')
            variant_name = (record.chrom[3:] if record.chrom.startswith('chr') else record.chrom, record.pos, record.ref, record.alts[0])
            variant = variants.get(variant_name, None)
            if variant is None:
                continue
            variant_csqs = set()
            lof = False
            for csq in record.info['CSQ']:
                csq = dict(zip(csq_header, csq.split('|')))
                if csq['BIOTYPE'] != 'protein_coding':
                    continue
                if csq['LoF'] == 'HC':
                    lof = True
                csqs = csq['Consequence'].split('&')
                if any(x in csqs for x in cds_variant_types):
                    variant_csqs.update(csqs)
            if not variant_csqs:
                variants.pop(variant_name)
                continue
            variant_most_severe_csq = None
            if 'splice_acceptor_variant' in variant_csqs or 'splice_donor_variant' in variant_csqs:
                variant_most_severe_csq = 'splice'
            elif 'stop_gained' in variant_csqs:
                variant_most_severe_csq = 'stop_gained'
            elif 'stop_lost' in variant_csqs:
                variant_most_severe_csq = 'stop_lost'
            elif 'start_lost' in variant_csqs:
                variant_most_severe_csq = 'start_lost'
            elif 'frameshift_variant' in variant_csqs:
                variant_most_severe_csq = 'frameshift'
            elif 'inframe_insertion' in variant_csqs:
                variant_most_severe_csq = 'inframe_insertion'
            elif 'inframe_deletion' in variant_csqs:
                variant_most_severe_csq = 'inframe_deletion'
            elif 'missense_variant' in variant_csqs:
                variant_most_severe_csq = 'missense'
            elif 'synonymous_variant' in variant_csqs or 'stop_retained_variant' in variant_csqs or 'start_retained_variant' in variant_csqs:
                variant_most_severe_csq = 'synonymous'
            else:
                print(f'WARNING (not in CDS): Variant {variant_name} ({variant_csqs}) will be omitted.')
                variants.pop(variant_name)
                continue
            categories = ['ALL', variant_most_severe_csq]
            if lof:
                categories.append('LOF')
            length = len(record.ref) - len(record.alts[0])
            if length == 0:
                if len(record.ref) > 1:
                    categories.append('MNP')
                else:
                    categories.append('SNP')
            elif length > 0:
                categories.append('INDEL')
                categories.append('DEL')
                length = abs(length)
                if length < 4:
                    categories.append(f'DEL:{length}')
                elif length < 10:
                    categories.append('DEL:4-9')
                else:
                    categories.append('DEL:10+')
            elif length < 0:
                categories.append('INDEL')
                categories.append('INS')
                length = abs(length)
                if length < 4:
                    categories.append(f'INS:{length}')
                elif length < 10:
                    categories.append('INS:4-9')
                else:
                    categories.append('INS:10+')
            variant.ac = record.info['AC'][0]
            variant.an = record.info['AN']
            variant.cat = categories


def load_dp_chunk(ifile, sample, chunk, chrom, pos):
    if chrom != chunk['chr'] or pos <= chunk['start'] or pos >= chunk['stop']:
        chunk['chr'] = chrom
        chunk['start'] = max(0, pos - 2)
        chunk['stop'] = pos + 500
        chunk['records'] = { r.pos: r.samples[sample]['DP'] for r in ifile.fetch(chunk['chr'], chunk['start'], chunk['stop']) }
        if not chunk['records']: # maybe we need to add 'chr' prefix
            chunk['records'] = { r.pos: r.samples[sample]['DP'] for r in ifile.fetch('chr' + chunk['chr'], chunk['start'], chunk['stop']) }



def read_depth(in_vcf1, in_vcf2, sample1, sample2, variants):
    with pysam.VariantFile(in_vcf1) as ifile1, pysam.VariantFile(in_vcf2) as ifile2:
        if sample1 not in set(ifile1.header.samples):
            raise Exception(f'Sample {sample1} was not found in {in_vcf1} file.')
        if sample2 not in set(ifile2.header.samples):
            raise Exception(f'Sample {sample2} was not found in {in_vcf2} file.')
        ifile1.subset_samples([sample1])
        ifile2.subset_samples([sample2])
        dp_chunk1 = { 'chr': '', 'start': -1, 'stop': -1, 'records': {} }
        dp_chunk2 = { 'chr': '', 'start': -1, 'stop': -1, 'records': {} }
        names_to_remove = []
        for name, (variant1, variant2) in variants.items(): # variants are in OrderedDict sorted by chrom, pos.
            chrom = name[0]
            pos = name[1]
            load_dp_chunk(ifile1, sample1, dp_chunk1, chrom, pos)
            load_dp_chunk(ifile2, sample2, dp_chunk2, chrom, pos)
            dp1 = dp_chunk1['records'].get(pos, None)
            dp2 = dp_chunk2['records'].get(pos, None)
            if (variant1.gt is not None and dp1 is None) or (variant2.gt is not None and dp2 is None):
                # if there was no entries in depth files, then most probably this region wasn't in our CDS list and we remove this variant completely from our comparison
                print(f'WARNING (no DP): Variant {name} will be omitted.')
                names_to_remove.append(name)
            else:
                variant1.dp = 0 if dp1 is None else dp1
                variant2.dp = 0 if dp2 is None else dp2
        for name in names_to_remove:
            variants.pop(name)


def get_variant_dp(in_vcf, sample, variant_name):
    #print('querying additional dp for', variant_name)
    with pysam.VariantFile(in_vcf) as ifile:
        vcf_samples = set((ifile.header.samples))
        if not sample in vcf_samples:
            raise Exception(f'Sample {sample} was not found in {in_vcf} file.')
        ifile.subset_samples([sample])
        chrom, pos, ref, alt = variant_name
        records = list(ifile.fetch(chrom, pos - 2, pos + 2)) # we fetch several positions and then filter out
        if not records: # maybe we need to add 'chr' prefix
           records = list(ifile.fetch('chr' + chrom, pos - 2, pos + 2))
        if not records: # at this point we assume that it was not sequenced
           return 0
        for record in records:
           if record.pos == pos:
              return record.samples[sample]['DP']
        return 0


def variant2summary(summary, variant):
    for category in variant.cat:
        summary['Count'].setdefault(category, 0)
        summary['Count'][category] += 1
        summary['DP'].setdefault(category, []).append(variant.dp)
        summary['AC'].setdefault(category, []).append(variant.ac)
        summary['AN'].setdefault(category, []).append(variant.an)
    if variant.dp > 0:
        for category in variant.cat:
            summary['Count_DP_0'].setdefault(category, 0)
            summary['Count_DP_0'][category] += 1
    if variant.dp > 20:
        for category in variant.cat:
            summary['Count_DP_20'].setdefault(category, 0)
            summary['Count_DP_20'][category] += 1


def summarize_sample(study, variants, qc_pass = True):
    summary = {
            'Count' : { },
            'Count_DP_0': { },
            'Count_DP_20': { },
            'DP': { },
            'AC': { },
            'AN': { }
            }
    for name, variant in ((x, y[study]) for x, y in variants.items()):
        if variant.gt is None:
            continue
        if qc_pass is not None and variant.is_pass != qc_pass:
            continue
        variant2summary(summary, variant)
    return summary


def summarize_pair(variants):
    summary = {}
    for name, (variant1, variant2) in variants.items():
        if variant1.gt is not None and variant2.gt is None:
            if not variant1.is_pass: # ignore if FAIL
                continue
            group = f'{variant1.is_pass}_{variant1.gt}_NONE_NONE'
            summary.setdefault(group, { 'study_1': { 'DP': {}, 'AN': {}, 'AC': {}}, 'study_2': { 'DP': {}}})
            for category in variant1.cat:
                summary[group]['study_1']['DP'].setdefault(category, []).append(variant1.dp)
                summary[group]['study_1']['AN'].setdefault(category, []).append(variant1.an)
                summary[group]['study_1']['AC'].setdefault(category, []).append(variant1.ac)
                summary[group]['study_2']['DP'].setdefault(category, []).append(variant2.dp)
        elif variant1.gt is None and variant2.gt is not None:
            if not variant2.is_pass: # ignore if FAIL
                continue
            group = f'NONE_NONE_{variant2.is_pass}_{variant2.gt}'
            summary.setdefault(group, { 'study_1': {'DP': {}}, 'study_2': { 'DP': {}, 'AN': {}, 'AC': {}}})
            for category in variant2.cat:
                summary[group]['study_1']['DP'].setdefault(category, []).append(variant1.dp)
                summary[group]['study_2']['DP'].setdefault(category, []).append(variant2.dp)
                summary[group]['study_2']['AN'].setdefault(category, []).append(variant2.an)
                summary[group]['study_2']['AC'].setdefault(category, []).append(variant2.ac)
        else:
            # both variants were found
            if not variant1.is_pass and not variant2.is_pass: # ignore if both FAIL
                continue
            group = f'{variant1.is_pass}_{variant1.gt}_{variant2.is_pass}_{variant2.gt}'
            summary.setdefault(group, { 'study_1': { 'DP': {}, 'AN': {}, 'AC': {}}, 'study_2': { 'DP': {}, 'AN': {}, 'AC': {}}})
            for category in variant1.cat: # since ref and alt alleles match, then categories should be the same in both studies
                summary[group]['study_1']['DP'].setdefault(category, []).append(variant1.dp)
                summary[group]['study_1']['AN'].setdefault(category, []).append(variant1.an)
                summary[group]['study_1']['AC'].setdefault(category, []).append(variant1.ac)
                summary[group]['study_2']['DP'].setdefault(category, []).append(variant2.dp)
                summary[group]['study_2']['AN'].setdefault(category, []).append(variant2.an)
                summary[group]['study_2']['AC'].setdefault(category, []).append(variant2.ac)
    return summary


def summarize(args, variants_union, output_suffix = None):
    with gzip.open(args.out_file + ('.' + output_suffix if output_suffix else '') +  '.json.gz', 'wt') as ofile:
        summaries = []
        for filter_value in [None, True, False]:
            summary = summarize_sample(0, variants_union, filter_value)
            summary['sample'] = args.sample1
            summary['gt_file'] = args.in_vcf1
            summary['dp_file'] = args.in_dp_vcf1
            summary['vep_file'] = args.in_vep_vcf1
            summary['study'] = 1
            summary['filter'] = filter_value
            summaries.append(summary)
        for filter_value in [None, True, False]:
            summary = summarize_sample(1, variants_union, filter_value)
            summary['sample'] = args.sample1
            summary['gt_file'] = args.in_vcf2
            summary['dp_file'] = args.in_dp_vcf2
            summary['vep_file'] = args.in_vep_vcf2
            summary['study'] = 2
            summary['filter'] = filter_value
            summaries.append(summary)
        summary = summarize_pair(variants_union)
        summary['pairwise'] = 'study1_vs_study2'
        summaries.append(summary)
        json.dump(summaries, ofile, sort_keys = True)


if __name__ == '__main__':
    args = argparser.parse_args()

    targets = dict()
    read_targets(args.in_targets_bed, targets)
    print(f'# Target regions:')
    for chrom, chrom_targets in targets.items():
        print(f'- {chrom}: {len(chrom_targets)} intervals over {sum(interval.end - interval.begin for interval in chrom_targets)} bp')

    variants1 = dict()
    t_start = time.time()
    read_genotypes(args.in_vcf1, args.sample1, variants1)
    print(f'# All variants in individual from study 1: {len(variants1)} ({time.time() - t_start} sec)')
    t_start = time.time()
    read_vep(args.in_vep_vcf1, variants1)
    print(f'# All coding variants in individual from study 1: {len(variants1)} ({time.time() - t_start} sec)')

    variants2 = dict()
    t_start = time.time()
    read_genotypes(args.in_vcf2, args.sample2, variants2)
    print(f'# All variants in individual from study 2: {len(variants2)} ({time.time() - t_start} sec)')
    t_start = time.time()
    read_vep(args.in_vep_vcf2, variants2)
    print(f'# All coding variants in individual from study 2: {len(variants2)} ({time.time() - t_start} sec)')

    t_start = time.time()
    variants_union = OrderedDict()
    for name in sorted(variants1.keys() | variants2.keys()):
        variants_union[name] = (
            variants1.pop(name, Variant(gt = None, is_pass = None)),
            variants2.pop(name, Variant(gt = None, is_pass = None))
        )
    print(f'# Union: {len(variants_union)} ({time.time() - t_start} sec)')

    t_start = time.time()
    read_depth(args.in_dp_vcf1, args.in_dp_vcf2, args.sample1, args.sample2, variants_union)
    print(f'# Union with DP: {len(variants_union)} ({time.time() - t_start} sec)')

    t_start = time.time()
    summarize(args, variants_union)
    print(f'# Summrized in {time.time() - t_start} sec.')

    t_start = time.time()
    target_variants_union = OrderedDict()
    while variants_union:
        name, (variant1, variant2) = variants_union.popitem(last = False)
        chrom = name[0]
        start = name[1]
        stop = start + max(len(name[2]), len(name[3])) # we don't substract 1 here because IntervalTree assumes open ended interval
        if targets[chrom].overlaps(start, stop):
            target_variants_union[name] = (variant1, variant2)
    print(f'# Union with DP on target: {len(target_variants_union)} ({time.time() - t_start} sec)')

    t_start = time.time()
    summarize(args, target_variants_union, 'on_target')
    print(f'# Summrized in {time.time() - t_start} sec.')

    print('Done.')
