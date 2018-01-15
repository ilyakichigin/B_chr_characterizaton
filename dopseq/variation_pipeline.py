#!/usr/bin/env python

import os
import sys
import yaml
import pysam
import argparse

from dopseq.tools import utils, \
                         vca_stat


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=
                    """Variant calling, filtering and annotation for
                    specified genome regions.

                    See config for process description, inputs and outputs.
                    """
                    )
    parser.add_argument('config_file', help='input configuration file')
    parser.add_argument('-c', '--create_makefile',
                        action='store_true', default=False,
                        help='create example makefile')
    parser.add_argument('-d', '--dry_run', action='store_true', default=False,
                        help='Check all dependencies and print all commands')

    return parser.parse_args()


def add_read_groups(bam, gatk, samtools, dry_run=True):
    """Add read group to BAM if not present"""

    rg_present = False

    with pysam.AlignmentFile(bam) as in_bam:
        if 'RG' in in_bam.header.keys():
            rg_present = True
    if not rg_present:
        bam_prefix = bam.split('.')[0]
        temp_bam = bam + '.rg.tmp'
        rg_args = ('%s AddOrReplaceReadGroups'
                   ' --RGLB 1 --RGPL Illumina --RGPU 1 --RGSM %s'
                   ' -I %s -O %s' % (gatk, bam_prefix, bam, temp_bam))
        utils.run_command(rg_args, dry_run=dry_run)
        sys.stderr.write('mv %s %s\n' % (temp_bam, bam))
        if not dry_run:
            os.rename(temp_bam, bam)


def generate_names(sample, ref, reg):
    """Output filenames generator"""
    genome = ref.split('/')[-1].split('.')[0]
    base = '%s/%s.%s' % (sample, sample, genome)
    assert reg.endswith('.bed')
    var_base = sample + '/' + reg.split('/')[-1][:-4]

    names = {
        # input files
        'bam': base + '.filter.bam',
        'pos': base + '.pos.bed',
        # output files
        'log': var_base + '.log',
        'reg_pos': var_base + '.pos.bed',
        'fadict': '.'.join(ref.split('.')[:-1]) + '.dict',
        'faidx': ref + '.fai',
        'md_bam': base + '.filter.rmdup.bam',
        'md_log': base + '.filter.rmdup.log',
        'vcf': var_base + '.vcf',
        'ann_vcf': var_base + '.ann.vcf',
        'ann_html': var_base + '.ann.html',
        'ann_stat': var_base + '.ann.stat.txt'
        }
    return names


def main():

    args = parse_command_line_arguments()

    if args.create_makefile:
        examples_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                    'examples'))
        utils.copy_to_cwd(examples_dir,
                          'variation_makefile.yaml', args.config_file)
        sys.stderr.write('Copied example makefile for variation pipeline to '
                         '%s.\nYou can now fill it with your data and '
                         'parameters, then test ("-d" option) and run the '
                         'pipeline.\n' % (args.config_file))
        sys.exit(0)

    if args.dry_run:
        sys.stderr.write('This is a dry run. Only command listing '
                         'will be produced.\n')

    with open(args.config_file) as conf_file:
        conf = yaml.safe_load(conf_file)
    for f in (conf['ref_genome'], conf['snpEff_path']):
        utils.check_file(f)
    if 'bedtools_path' not in conf.keys():
        conf['bedtools_path'] = 'bedtools'
    for p in (conf['bedtools_path'], conf['gatk_path']):
        utils.check_exec(p)

    for sample in conf['samples']:

        sys.stderr.write('----Processing sample %s with variation pipeline'
                         '----\n' % sample)

        # generate filenames and check input files
        reg_bed = conf['samples'][sample]
        fnames = generate_names(sample, conf['ref_genome'], reg_bed)
        utils.check_file(reg_bed)
        utils.check_file(fnames['bam'])

        # create sequence dictionary for genome
        if not os.path.isfile(fnames['fadict']):
            csd_args = '%s CreateSequenceDictionary -R %s -O %s' % (
                    conf['gatk_path'],
                    conf['ref_genome'],
                    fnames['fadict'])
            utils.run_command(csd_args, errfile=fnames['log'],
                              dry_run=args.dry_run)

        # index genome
        if not os.path.isfile(fnames['faidx']):
            sys.stderr.write('Indexing %s with pysam\n' % conf['ref_genome'])
            if not args.dry_run:
                pysam.faidx(conf['ref_genome'])

        # create BED file for positions within regions
        if not os.path.isfile(fnames['reg_pos']):
            rp_args = '%s intersect -a %s -b %s' % (conf['bedtools_path'],
                                                    fnames['pos'],
                                                    reg_bed)
            utils.run_command(rp_args,
                              outfile=fnames['reg_pos'], errfile=fnames['log'],
                              dry_run=args.dry_run)
        else:
            sys.stderr.write('%s bed file with read positions within regions'
                             'exists. OK!\n'% fnames['reg_pos'])

        # mark duplicates in BAM file
        if not os.path.isfile(fnames['md_bam']):
            add_read_groups(fnames['bam'],
                            conf['gatk_path'],
                            conf['samtools_path'],
                            dry_run=args.dry_run)
            md_args = '%s MarkDuplicates -I %s -O %s -M %s' % (
                    conf['gatk_path'],
                    fnames['bam'],
                    fnames['md_bam'],
                    fnames['md_log'])
            utils.run_command(md_args, errfile=fnames['log'],
                              dry_run=args.dry_run)
        else:
            sys.stderr.write('%s deduplicated bam file index exists. OK!\n'
                             % fnames['md_bam'])

        # index dedulpicated BAM file
        if not os.path.isfile(fnames['md_bam'] + '.bai'):
            sys.stderr.write('Indexing %s with pysam\n' % fnames['md_bam'])
            if not args.dry_run:
                pysam.index(fnames['md_bam'])


        # infer genotypes for verified regions only
        if not os.path.isfile(fnames['vcf']):
            hc_args = '%s HaplotypeCaller -I %s -O %s -R %s -L %s' % (
                    conf['gatk_path'],
                    fnames['md_bam'],
                    fnames['vcf'],
                    conf['ref_genome'],
                    fnames['reg_pos'])
            utils.run_command(hc_args, errfile=fnames['log'], dry_run=args.dry_run)
        else:
            sys.stderr.write('%s file with region variants exists. OK!\n'
                             % fnames['vcf'])

        # annotate variants
        if not os.path.isfile(fnames['ann_vcf']):
            se_args = 'java -Xmx%s -jar %s -s %s %s %s' % (
                    conf['max_mem'],
                    conf['snpEff_path'],
                    fnames['ann_html'],
                    conf['snpEff_genome'],
                    fnames['vcf'])
            utils.run_command(se_args,
                              outfile=fnames['ann_vcf'], errfile=fnames['log'],
                              dry_run=args.dry_run)
        else:
            sys.stderr.write('%s file with annotated variants exists. OK!\n'
                             % fnames['ann_vcf'])

        # calculate annotation statistics
        if not os.path.isfile(fnames['ann_stat']):
            vs_args = argparse.Namespace(
                                         annotated_vcf=fnames['ann_vcf'],
                                         reg_bed=reg_bed,
                                         regpos_bed=fnames['reg_pos'],
                                         snpEff_path=conf['snpEff_path'],
                                         genome=conf['snpEff_genome'],
                                         max_mem=conf['max_mem'],
                                         bedtools_path=conf['bedtools_path'],
                                         out_file=fnames['ann_stat'],
                                         log_file=fnames['log'],
                                         dry_run=args.dry_run)
            vca_stat.main(vs_args)
        else:
            sys.stderr.write('%s file with annotation statistics exists. OK!\n'
                             % fnames['ann_stat'])

        sys.stderr.write('----Done processing sample %s with variation '
                         'pipeline----\n\n' % sample)

    if args.dry_run:
        sys.stderr.write('This was a dry run. Only command listing was '
                         'produced.\n')


if __name__ == '__main__':
    main()
