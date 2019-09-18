import argparse
import pysam
import statistics
import json


argparser = argparse.ArgumentParser(description = 'Compares called genotypes in two samples from different sequencing experiments')
argparser.add_argument('-s1', '--sample-study1', metavar = 'name', dest = 'sample1', required = True, help = 'Sample name from study 1.')
argparser.add_argument('-g1', '--genotypes-study1', metavar = 'file', dest = 'in_vcf1', required = True, help = 'Genotypes from study 1 in VCF/BCF.')
argparser.add_argument('-d1', '--depth-study1', metavar = 'file', dest = 'in_dp_vcf1', required = True, help = 'Depth (DP) from study 1 in VCF/BCF.')
argparser.add_argument('-a1', '--annotations-study1', metavar = 'file', dest = 'in_vep_vcf1', required = True, help = 'VEP annotations from study 1 in VCF/BCF.')
argparser.add_argument('-s2', '--sample-study2', metavar = 'name', dest = 'sample2', required = True, help = 'Sample name from study 2.')
argparser.add_argument('-g2', '--genotypes-study2', metavar = 'file', dest = 'in_vcf2', required = True, help = 'Genotypes from study 2 in VCF/BCF.')
argparser.add_argument('-d2', '--depth-study2', metavar = 'file', dest = 'in_dp_vcf2', required = True, help = 'Depth (DP) from study 2 in VCF/BCF.')
argparser.add_argument('-a2', '--annotations-study2', metavar = 'file', dest = 'in_vep_vcf2', required = True, help = 'VEP annotations from study 2 in VCF/BCF.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_file', required = True, help = 'Output file name.')


cds_variant_types = ['start_lost', 'start_retained_variant', 'stop_retained_variant', 'synonymous_variant', 'missense_variant', 'inframe_insertion', 'inframe_deletion', 'stop_gained', 'stop_lost', 'frameshift_variant', 'coding_sequence_variant', 'splice_donor_variant', 'splice_acceptor_variant']


def read_vep(in_vep_vcf, variants):
    with pysam.VariantFile(in_vep_vcf) as ifile:
        csq_meta = ifile.header.info.get('CSQ', None)
        if csq_meta is None:
            raise Exception('No meta-information entry about CSQ INFO field found!')
        csq_header = csq_meta.description.split(':', 1)[1].strip().split('|')
        for record in ifile.fetch():
            variant_name = f'{record.chrom[3:] if record.chrom.startswith("chr") else record.chrom}-{record.pos}' # we work with bi-allelic variants only -> chr:pos is enough to identify variant
            variant_csqs = set()
            for csq in record.info['CSQ']:
                csq = dict(zip(csq_header, csq.split('|')))
                if csq['BIOTYPE'] != 'protein_coding':
                    continue
                csqs = csq['Consequence'].split('&')
                if any(x in csqs for x in cds_variant_types):
                    variant_csqs.update(csqs)
            if not variant_csqs:
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
                print(variant_name, variant_csqs)
                continue
            variants[variant_name] = { 'ref': record.ref, 'alt': record.alts[0], 'type': variant_most_severe_csq }


def read_genotypes(in_vcf, sample, variants):
    with pysam.VariantFile(in_vcf) as ifile:
        vcf_samples = set((ifile.header.samples))
        if not sample in vcf_samples:
            raise Exception(f'Sample {sample} was not found in {in_vcf} file.')
        ifile.subset_samples([sample])
        for record in ifile:
            variant_name = f'{record.chrom[3:] if record.chrom.startswith("chr") else record.chrom}-{record.pos}'
            if variant_name not in variants:
                continue
            gt = record.samples[sample]['GT']
            if None in gt or sum(gt) == 0:
                variants.pop(variant_name) # remove if alleles are missing or only reference allele was called
            else:
                variants[variant_name].update({ 'gt': sum(gt), 'pass': 'PASS' in list(record.filter)})


def read_depth(in_vcf, sample, variants):
    with pysam.VariantFile(in_vcf) as ifile:
        vcf_samples = set((ifile.header.samples))
        if not sample in vcf_samples:
            raise Exception(f'Sample {sample} was not found in {in_vcf} file.')
        ifile.subset_samples([sample])
        variants_to_remove = []
        for variant_name, variant_info in variants.items():
            chrom, pos = variant_name.split('-')
            pos = int(pos)
            records = list(ifile.fetch(chrom, pos - 2, pos + 2)) # we fetch several positions and then filter out
            if not records: # maybe we need to add 'chr' prefix
                records = list(ifile.fetch('chr' + chrom, pos - 2, pos + 2))
            if not records:
                # if there was no entries, then most probably this region wasn't in our CDS list and we remove this variant completely
                print(variant_name)
                variants_to_remove.append(variant_name)
                continue
            dp = None
            for record in records:
                if record.pos == pos:
                    dp = record.samples[sample]['DP']
                    break
            variant_info.update({ 'dp': 0 if dp is None else dp })
        for variant_name in variants_to_remove:
            variants.pop(variant_name)


def get_variant_categories(info):
    categories = ['ALL', info['type']]
    length = len(info['ref']) - len(info['alt'])
    if length == 0:
        if len(info['ref']) > 1:
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
    return categories


def summarize_sample(variants, qc_pass = True):
    summary = {
            'Count' : { },
            'Count_DP_0': { },
            'Count_DP_20': { },
            'DP': { }
            }
    for name, info in variants.items():
        if qc_pass is not None and info['pass'] != qc_pass:
            continue
        dp = info['dp']
        categories = get_variant_categories(info)
        for category in categories:
            summary['Count'].setdefault(category, 0)
            summary['Count'][category] += 1
            summary['DP'].setdefault(category, []).append(dp)
        if dp > 0:
            for category in categories:
                summary['Count_DP_0'].setdefault(category, 0)
                summary['Count_DP_0'][category] += 1
        if dp > 20:
            for category in categories:
                summary['Count_DP_20'].setdefault(category, 0)
                summary['Count_DP_20'][category] += 1
    return summary


def summarize_pair(variants1, variants2):
    names_union = set()
    names_union.update(variants1.keys())
    names_union.update(variants2.keys())

    summary = {}

    print('Union = ', len(names_union))
    for name in names_union:
        variant1 = variants1.get(name, None)
        variant2 = variants2.get(name, None)
        if variant1 is not None and variant2 is None:
            if not variant1['pass']: # ignore if FAIL
                continue
            group = f'{variant1["pass"]}_{variant1["gt"]}_NONE_NONE'
            summary.setdefault(group, { 'study_1': { 'DP': [] }})
            for category in get_variant_categories(variant1):
                summary[group]['study_1']['DP'].append(variant1['dp'])
        elif variant1 is None and variant2 is not None:
            if not variant2['pass']: # ignore if FAIL
                continue
            group = f'NONE_NONE_{variant2["pass"]}_{variant2["gt"]}'
            summary.setdefault(group, { 'study_2': {'DP': [] }})
            for category in get_variant_categories(variant2):
                summary[group]['study_2']['DP'].append(variant2['dp'])
        else:
            # both variants were found
            if not variant1['pass'] and not variant2['pass']: # ignore if both FAIL
                continue
            if variant1['ref'] == variant2['ref'] and variant1['alt'] == variant2['alt']: # both variants have same alllese
                group = f'{variant1["pass"]}_{variant1["gt"]}_{variant2["pass"]}_{variant2["gt"]}'
                summary.setdefault(group, { 'study_1': { 'DP': [] }, 'study_2': { 'DP': [] }})
                for category in get_variant_categories(variant1): # since ref and alt alleles match, then categories should be the same in both studies
                    summary[group]['study_1']['DP'].append(variant1['dp'])
                    summary[group]['study_2']['DP'].append(variant2['dp'])
            else:
                group = f'ALLELE_MISMATCH'
                summary.setdefault(group, { 'study_1': {}, 'study_2': {} })
                for category in get_variant_categories(variant1):
                    summary[group]['study_1'].setdefault(category, {'DP': []})
                    summary[group]['study_1'][category]['DP'].append(variant1['dp'])
                for category in get_variant_categories(variant2):
                    summary[group]['study_2'].setdefault(category, {'DP': []})
                    summary[group]['study_2'][category]['DP'].append(variant2['dp'])
    return summary


if __name__ == '__main__':
    args = argparser.parse_args()
    variants1 = dict()
    read_vep(args.in_vep_vcf1, variants1)
    print('All 1 = ', len(variants1))
    read_genotypes(args.in_vcf1, args.sample1, variants1)
    print('Sample 1 = ', len(variants1))
    read_depth(args.in_dp_vcf1, args.sample1, variants1)
    print('With DP 1 = ', len(variants1))

    variants2 = dict()
    read_vep(args.in_vep_vcf2, variants2)
    print('All 2 = ', len(variants2))
    read_genotypes(args.in_vcf2, args.sample2, variants2)
    print('Sample 2 = ', len(variants2))
    read_depth(args.in_dp_vcf2, args.sample2, variants2)
    print('With DP 2 = ', len(variants2))

    with open(args.out_file, 'w') as ofile:
        summaries = []
        for filter_value in [None, True, False]:
            summary = summarize_sample(variants1, filter_value)
            summary['sample'] = args.sample1
            summary['study'] = 1
            summary['filter'] = filter_value
            summaries.append(summary)
        for filter_value in [None, True, False]:
            summary = summarize_sample(variants2, filter_value)
            summary['sample'] = args.sample1
            summary['study'] = 2
            summary['filter'] = filter_value
            summaries.append(summary)
        summary = summarize_pair(variants1, variants2)
        summary['pairwise'] = 'study1_vs_study2'
        summaries.append(summary)
        json.dump(summaries, ofile, indent = 4)

