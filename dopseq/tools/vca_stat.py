#!/usr/bin/env python

import os
import sys
import argparse

from collections import defaultdict, OrderedDict

from dopseq.tools import utils


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=
                    """
                    Calculate variant annotation statistics: top variant
                    effects taken from annotated_vcf, sequence revovery -
                    from regpos_bed, gene list - from reg_bed.
                    """
                    )
    parser.add_argument('annotated_vcf', help='vcf file annotated with snpEff')
    parser.add_argument('reg_bed', help='bed file with target regions')
    parser.add_argument('regpos_bed', help='bed file with positions covered '
                                            'by reads within target regions')
    parser.add_argument('-s', '--snpEff_path', help='path to snpEff jar file')
    parser.add_argument('-g', '--genome', help='snpEff genome version used '
                                               'for annotation')
    parser.add_argument('-m', '--max_mem', default='4g',
                        help='maximum RAM allocated to snpEff. Default: 4 Gb')
    parser.add_argument('-b', '--bedtools_path', help='path to bedtools')
    parser.add_argument('-o', '--out_file', help='output statistics file')
    parser.add_argument('-l', '--log_file', help='log file')
    parser.add_argument("-d", "--dry_run", default=False,
                        help="Perform only a dry run")

    return parser.parse_args()


def generate_filenames(bed):
    """Generate filenames for temporary intermediate files"""

    assert bed.endswith('.bed')
    bed_base = bed.split('/')[-1][:-4]

    names = {
        'regpos_features_txt': 'temp.%s.pos.features.txt' % bed_base,
        'regpos_features_summary_txt': ('temp.%s.pos.features.summary.txt'
                                        % bed_base),
        'regpos_features_summary_html': ('temp.%s.pos.features.summary.html'
                                         % bed_base),
        'reg_features_txt': 'temp.%s.features.txt' % bed_base,
        'reg_features_summary_txt': ('temp.%s.features.summary.txt'
                                     % bed_base),
        'reg_features_summary_html': ('temp.%s.features.summary.html'
                                      % bed_base),
        'intron_bed': 'temp.%s.intron.bed' % bed_base,
        'one_reg': 'temp.%s.onereg.bed' % bed_base
        }

    return names


def variant_types():
    """Conevsion dict of variant types from snpEff to simple categories:
    'cds', 'exon' (except for cds), 'intron', 'intergenic', and 'other' -
    for unclassified variants. Variant types taken from 
    https://github.com/pcingola/SnpEff/blob/919c61cef2397ef660db63a3ffd1c902f0cf3a40/src/main/java/org/snpeff/snpEffect/EffectType.java
    """

    return {'coding_sequence_variant': 'cds',
            'inframe_insertion': 'cds',
            'conservative_inframe_insertion': 'cds',
            'disruptive_inframe_insertion': 'cds',
            'inframe_deletion': 'cds',
            'conservative_inframe_deletion': 'cds',
            'disruptive_inframe_deletion': 'cds',
            'frameshift_variant': 'cds',
            'missense_variant': 'cds',
            'initiator_codon_variant': 'cds',
            'non_canonical_start_codon': 'cds',
            'stop_retained_variant': 'cds',
            'protein_protein_contact': 'cds',
            'structural_interaction_variant': 'cds',
            'rare_amino_acid_variant': 'cds',
            'stop_lost': 'cds',
            'start_lost': 'cds',
            'stop_gained': 'cds',
            'synonymous_variant': 'cds',
            'start_retained': 'cds',

            'exon_region': 'exon',
            'exon_variant': 'exon',
            'exon_loss_variant': 'exon',
            'miRNA': 'exon',
            '5_prime_UTR_premature_start_codon_gain_variant': 'exon',
            '5_prime_UTR_premature_start_codon_variant': 'exon',
            '5_prime_UTR_variant': 'exon',
            '5_prime_UTR_truncation': 'exon',
            '3_prime_UTR_variant': 'exon',
            '3_prime_UTR_truncation': 'exon',
            'non_coding_transcript_exon_variant': 'exon',
            'non_coding_transcript_variant': 'exon',

            'intron_variant': 'intron',
            'conserved_intron_variant': 'intron',
            'splice_acceptor_variant': 'intron',
            'splice_donor_variant': 'intron',
            'splice_region_variant': 'intron',
            'splice_branch_variant': 'intron',

            'downstream_gene_variant': 'intergenic',
            'intergenic_region': 'intergenic',
            'conserved_intergenic_variant': 'intergenic',
            'regulatory_region_variant': 'intergenic',
            'upstream_gene_variant': 'intergenic',
            'TF_binding_site_variant': 'intergenic',
            'TFBS_ablation': 'intergenic',

            'chromosome': 'other',
            'duplication': 'other',
            'inversion': 'other',
            'exon_loss': 'other',
            'exon_loss_variant': 'other',
            'exon_region': 'other',
            'gene_variant': 'other',
            'feature_ablation': 'other',
            'feature_elongation': 'other',
            'feature_fusion': 'other',
            'gene_fusion': 'other',
            'bidirectional_gene_fusion': 'other',
            'rearranged_at_DNA_level': 'other',
            'intragenic_variant': 'other',
            'transcript_variant': 'other',
            'transcript_ablation': 'other',
            'rearranged_at_DNA_level': 'other',
            'sequence_feature': 'other',
            'custom': 'other'
            }


def init_tr_features(gene_id, type_id):
    '''Return dict with transcript features used in analysis'''

    return {'Gene_name': '', # from VCF - unneeded
            'Gene_id': gene_id,
            'Transcript_type': type_id,
            'Transcript_coord': [], # ['chr\tstart\tend\n']
            'Exon_coord': [],
            'CDS_coord': [],
            'Intron_coord': [],
            'Transcript_span': [], # [bp]
            'CDS_bps': [], # [bp]
            'Exon_bps': defaultdict(list), # exon_no: [bp]
            'Intron_bps': defaultdict(list), # intron_no: [bp]
            'CDS_bp_sum': 0,
            'Exon_bp_sum': 0,
            'Intron_bp_sum': 0,
            'Transcript_bp_sum': 0
            }


def var_eff_counts(handle):
    """Return counts of top variant effects (cds>exon>intron>intergenic)
    for file handle or list of lines"""

    eff_count = {'cds': 0, 'exon': 0, 'intron': 0, 'intergenic': 0, 'other': 0}
    unrecognized_types = []
    v_types = variant_types()
    for l in handle:
        if l.startswith('#'):
            continue
        ll = l.split()
        if len(ll) < 10:
            continue
        v_info = ll[7].split('|')
        v_effects = []
        # iterate over multiple effects
        for i in range(len(v_info) / 15):
            for t in v_info[i * 15 + 1].split('&'):
                try:
                    eff = v_types[t]
                # variant type not in variant_types()
                except:
                    unrecognized_types.append(t)
                    eff = 'other'
            v_effects.append(eff)
        if 'cds' in v_effects:
            eff_count['cds'] += 1
        elif 'exon' in v_effects:
            eff_count['exon'] += 1
        elif 'intron' in v_effects:
            eff_count['intron'] += 1
        elif 'intergenic' in v_effects:
            eff_count['intergenic'] += 1
        else:
            eff_count['other'] += 1
    if len(unrecognized_types) > 0:
        sys.stderr.write('Unrecognized variant types:\n'
                         + '\n'.join(set(unrecognized_types)) + '\n')

    return eff_count


def bed_size(handle):
    """Return total size of regions in bed file handle or list of lines"""
    
    size = 0
    for l in handle:
        ll = l.split()
        if len(ll) < 3:
            continue
        size += int(ll[2]) - int(ll[1])
    
    return size


def main(args):

    fnames = generate_filenames(args.reg_bed)

    # pre-clean up
    if not args.dry_run:
        for f in fnames:
            if os.path.isfile(f):
                os.unlink(fnames[f])

    # generate features file for regpos and reg intervals
    for (in_bed, out_cnt) in ((args.regpos_bed, fnames['regpos_features_txt']),
                              (args.reg_bed, fnames['reg_features_txt'])):
        if not os.path.isfile(out_cnt):
            feature_args = 'java -Xmx%s -jar %s count -n %s %s %s' % (
                    args.max_mem,
                    args.snpEff_path,
                    out_cnt[:-4],
                    args.genome,
                    in_bed)
            utils.run_command(feature_args, errfile=args.log_file,
                              dry_run=args.dry_run)

    if args.dry_run:
        return 0

    # bp recovery data from features file
    reg_fs = dict()
    for features_file in (fnames['regpos_features_txt'],
                          fnames['reg_features_txt']):
        extended = (True if features_file == fnames['regpos_features_txt']
                    else False)
        with open(features_file) as f:
            for l in f:
                if l.startswith('chr\tstart'):
                    continue
                ll = l.split()
                if len(ll) < 6:
                    continue
                # convert to 0-based coordiates for BED
                ll[1] = str(int(ll[1]) - 1)
                coordinates = '\t'.join(ll[0:3]) + '\n'
                span = int(ll[2]) - int(ll[1])
                feature_data = ll[3].split(';')
                feature_data = [data.split(':') for data in feature_data]
                bp = int(ll[5])
                feature_type = feature_data[0][0]
                if feature_type == 'Transcript':
                    transcript = feature_data[0][1]
                    if transcript not in reg_fs.keys():
                        reg_fs[transcript] = init_tr_features(
                                feature_data[1][1],
                                feature_data[0][2])
                    # regpos transcripts interpreted as CDS
                    if extended:
                        reg_fs[transcript]['CDS_bps'].append(bp)
                        reg_fs[transcript]['CDS_coord'].append(coordinates)
                    # reg-only transcripts interpreted as transcripts
                    else:
                        reg_fs[transcript]['Transcript_coord'].append(coordinates)
                    reg_fs[transcript]['Transcript_span'].append(span)
                elif feature_type in ('Exon','Intron') and extended:
                    transcript = feature_data[1][1]
                    if transcript not in reg_fs.keys():
                        reg_fs[transcript] = init_tr_features(
                                feature_data[2][1],
                                feature_data[1][2])
                    feature_no = feature_data[0][1]
                    reg_fs[transcript][feature_type + '_bps'][feature_no].append(bp)
                    reg_fs[transcript][feature_type + '_coord'].append(coordinates)

    # per transcript bp recovery data
    for transcript in reg_fs.keys():

        transcript_span = max(reg_fs[transcript]['Transcript_span'])
        i = reg_fs[transcript]['Transcript_span'].index(transcript_span)
        # transcript listed in regpos
        if len(reg_fs[transcript]['CDS_coord']) > 0:
            reg_fs[transcript]['Transcript_coord'] = ('%s\t%s\n' %
                    (reg_fs[transcript]['CDS_coord'][i][:-1], transcript))
            reg_fs[transcript]['Transcript_bp_sum'] = reg_fs[transcript]['CDS_bps'][i]
            del reg_fs[transcript]['CDS_coord'][i]
            del reg_fs[transcript]['CDS_bps'][i]
            cds_bp = sum(reg_fs[transcript]['CDS_bps'])
            reg_fs[transcript]['CDS_bp_sum'] = cds_bp
            
            exon_bp = 0
            for exon in reg_fs[transcript]['Exon_bps'].keys():
                exon_bp += max(reg_fs[transcript]['Exon_bps'][exon])
            reg_fs[transcript]['Exon_bp_sum'] = exon_bp

            intron_bp = 0
            for intron in reg_fs[transcript]['Intron_bps'].keys():
                intron_bp += max(reg_fs[transcript]['Intron_bps'][intron])
            reg_fs[transcript]['Intron_bp_sum'] = intron_bp

            # debug
            if reg_fs[transcript]['Transcript_bp_sum'] != exon_bp + intron_bp:
                print transcript
                print reg_fs[transcript]
                sys.exit()
        # transcript in reg, but not in regpos
        else:
            reg_fs[transcript]['Transcript_coord'] = ('%s\t%s\n' %
                    (reg_fs[transcript]['Transcript_coord'][i][:-1], transcript))
            

    # beds with feature recovery
    beds = {
            'merged_cds': '',
            'merged_exon': '',
            'merged_intron': '',
            'shrinked_intron': '',
            'all_transcript': ''.join(reg_fs[tr]['Transcript_coord'] 
                                      for tr in reg_fs.keys())
            }
    for (feature, out_feature) in (('CDS_coord', 'merged_cds'),
                                  ('Exon_coord', 'merged_exon'),
                                  ('Intron_coord', 'merged_intron')):
        feature_coord = ''
        for transcript in reg_fs.keys():
            if 'Transcript_bp_sum' > 0:
                feature_coord += ''.join(reg_fs[transcript][feature])
        sorted_bed = utils.run_command(
                args.bedtools_path + ' sort -i -',
                stdin=feature_coord,
                verbose=False, dry_run=False,
                return_out=True)
        merged_bed = utils.run_command(
                args.bedtools_path + ' merge -i -',
                stdin=sorted_bed,
                verbose=False, dry_run=False,
                return_out=True)
        merged_regpos_bed = utils.run_command(
                '%s intersect -a - -b %s' % (args.bedtools_path,
                                            args.regpos_bed),
                stdin=merged_bed,
                verbose=False, dry_run=False,
                return_out=True)
        beds[out_feature] = merged_regpos_bed
    with open(fnames['intron_bed'],'w') as f:
        f.write(beds['merged_intron'])
    beds['shrinked_intron'] = utils.run_command(
                '%s subtract -a %s -b -' % (args.bedtools_path,
                                            fnames['intron_bed']),
                stdin=beds['merged_exon'],
                verbose=False, dry_run=False,
                return_out=True)

    # top variant effects from annotated vcf
    with open(args.annotated_vcf) as f:
        eff_count = var_eff_counts(f)

    # summary statistics
    stats = [['#feature','bp_covered','variants','variants_per_kbp'],
             ['total','','',''],
             ['cds','','',''],
             ['exons','','',''],
             ['introns','','',''],
             ['intergenic','','','']]
    # coverage statistics
    stats[1][1] = bed_size(open(args.regpos_bed))
    stats[2][1] = bed_size(beds['merged_cds'].split('\n'))
    stats[3][1] = bed_size(beds['merged_exon'].split('\n'))
    stats[4][1] = bed_size(beds['shrinked_intron'].split('\n'))
    stats[5][1] = stats[1][1] - stats[3][1] - stats[4][1]
    # variant statistics
    stats[1][2] = (eff_count['cds'] + eff_count['exon']
                  + eff_count['intron'] + eff_count['intergenic'])
    stats[2][2] = eff_count['cds']
    stats[3][2] = eff_count['exon'] + eff_count['cds']
    stats[4][2] = eff_count['intron']
    stats[5][2] = eff_count['intergenic']
    # variant per kbp statistics
    for i in range(1, len(stats)):
        try:
            stats[i][3] = '%.2f' % (float(stats[i][2]) / stats[i][1] * 1000)
        except:
            stats[i][3] = 'NA'
    for i in range(1, len(stats)):
        for j in range(1, len(stats[0]) - 1):
            stats[i][j] = str(stats[i][j])
    
    # per region statistics
    reg_stats = ('#region\tbp_total\tbp_exon\tbp_intron\tbp_intergenic\t'
                 'var_total\tvar_exon\tvar_intron\tvar_intergenic\t'
                 'var_per_kbp_total\tvar_per_kbp_exon\tvar_per_kbp_intron\t'
                 'var_per_kbp_intergenic\ttranscripts\n')
    with open(args.reg_bed) as f:
        for l in f:
            with open(fnames['one_reg'], 'w') as o:
                o.write(l)
            reg_eff_count = var_eff_counts(utils.run_command(
                    '%s intersect -a %s -b %s' % (args.bedtools_path,
                                                  args.annotated_vcf,
                                                  fnames['one_reg']),
                    verbose=False, dry_run=False, return_out=True).split('\n'))
            total_bp = bed_size(utils.run_command(
                    '%s intersect -a %s -b %s' % (args.bedtools_path,
                                                  args.regpos_bed,
                                                  fnames['one_reg']),
                    verbose=False, dry_run=False, return_out=True).split('\n'))
            exon_bp = bed_size(utils.run_command(
                    '%s intersect -a - -b %s' % (args.bedtools_path,
                                                 fnames['one_reg']),
                    stdin=beds['merged_exon'],
                    verbose=False, dry_run=False, return_out=True).split('\n'))
            intron_bp = bed_size(utils.run_command(
                    '%s intersect -a - -b %s' % (args.bedtools_path,
                                                 fnames['one_reg']),
                    stdin=beds['shrinked_intron'],
                    verbose=False, dry_run=False, return_out=True).split('\n'))
            intergenic_bp = total_bp - exon_bp - intron_bp
            reg_transcripts = utils.run_command(
                    '%s intersect -a - -b %s' % (args.bedtools_path,
                                                  fnames['one_reg']),
                    stdin=beds['all_transcript'],
                    verbose=False, dry_run=False, return_out=True).split('\n')
            reg_transcript_ids = []
            if len(reg_transcripts) > 0:
                for rt in reg_transcripts:
                    rtd = rt.split('\t')
                    if len(rtd) > 3:
                        reg_transcript_ids.append(rtd[3])
            
            ll = l.split()
            try:
                total_vpk = '%.2f' % (sum(reg_eff_count.values())
                                      /float(total_bp)*1000)
            except:
                total_vpk = 'NA'
            try:
                exon_vpk = '%.2f' % ((reg_eff_count['cds'] + reg_eff_count['exon'])
                                     /float(exon_bp)*1000)
            except:
                exon_vpk = 'NA'
            try:
                intron_vpk = '%.2f' % (reg_eff_count['intron']
                                       /float(intron_bp)*1000)
            except:
                intron_vpk = 'NA'
            try:
                intergenic_vpk = '%.2f' % (reg_eff_count['intergenic']
                                           /float(intergenic_bp)*1000)
            except:
                intergenic_vpk = 'NA'

            reg_stats += ('%s:%s-%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t'
                          '%s\t%s\t%s\t%s\t%s\n' % 
                    (ll[0], ll[1], ll[2], total_bp, exon_bp, intron_bp, 
                     intergenic_bp, sum(reg_eff_count.values()),
                     reg_eff_count['cds'] + reg_eff_count['exon'], 
                     reg_eff_count['intron'], reg_eff_count['intergenic'],
                     total_vpk, exon_vpk, intron_vpk, intergenic_vpk,
                     ';'.join(reg_transcript_ids)))

    # write summary file
    genes_file = args.annotated_vcf[:-3] + 'genes.txt'
    processed_transcripts = []
    novar_transcripts = []
    outside_transcripts = []
    with open(genes_file) as f, \
         open(args.out_file, 'w') as o:

        o.write('# Variants:\t%s\n# Regions:\t%s\n\n\n' % (args.annotated_vcf,
                                                           args.reg_bed))

        o.write('# Summary statistics\n')
        for s in stats:
            o.write('\t'.join(s) + '\n')
        o.write('\n%d unclassified variants not included\n\n\n'
                % eff_count['other'])

        o.write('# Per region statistics\n')
        o.write(reg_stats)
        o.write('\n\n')

        o.write('# Per gene statistics\n')
        # genes listed by snpEff
        for l in f:
            ll = l.split('\t')
            if l.startswith('# '):
                #o.write(l)
                continue
            elif l.startswith('#GeneName'):
                ll = l.split('\t')
                ll = ll[:4] + ['Exon_bp_sum', 'Intron_bp_sum'] + ll[4:]
                o.write('\t'.join(ll))
                continue
            elif len(ll) < 4:
                #o.write(l)
                continue
            try:
                processed_transcripts.append(ll[2])
                bp_data = [str(reg_fs[ll[2]]['Exon_bp_sum']),
                           str(reg_fs[ll[2]]['Intron_bp_sum'])]
                o.write('\t'.join(ll[:4] + bp_data + ll[4:]))
            except:
                outside_transcripts.append(ll[:4] + ['0', '0'] + ll[4:])

        # transcripts not listed by snpEff - no variants
        for transcript in reg_fs.keys():
            if transcript not in processed_transcripts:
                processed_transcripts.append(transcript)
                fdata = reg_fs[transcript]
                novar_transcripts.append('NA\t%s\t%s\t%s\t%d\t%d\t'
                        '0\t0\t0\t0' % (fdata['Gene_id'],
                                        transcript,
                                        fdata['Transcript_type'],
                                        fdata['Exon_bp_sum'],
                                        fdata['Intron_bp_sum']))
        novar_transcripts.sort()
        o.write('\n'.join(novar_transcripts) + '\n\n')

        o.write('# Genes outside regions\n')
        for tr in outside_transcripts:
            o.write('\t'.join(tr))

    # clean up
    for f in fnames:
        if os.path.isfile(f):
            os.unlink(fnames[f])    


if __name__ == '__main__':
    main(parse_command_line_arguments())