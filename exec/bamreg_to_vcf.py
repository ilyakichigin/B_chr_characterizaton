#!/usr/bin/env python

# From bam file create annotated VCF callset for target genomic regions and statistics. See .conf for details.
# Requires: GATK (tested on v.3.3.0), picard-tools (tested on v.1.125)
# Requires in $PATH: samtools (in $PATH, tested on v.0.1.19), bedtools
# Usage: bamreg_to_vcf.py bamreg_to_vcf.conf

import subprocess
import sys
import ConfigParser
import os

def prepare_genome(genome_fasta, path_to_picard):

    # Using picard tools prepare genomic fasta to be used in GATK    
    
    print '-----prepare-genome'
    dict_file = '.'.join(genome_fasta.split('.')[:-1])+'.dict'
    if not os.path.isfile(dict_file):
        pg_command = ['java','-Xmx1g','-jar',path_to_picard,
                      'CreateSequenceDictionary',
                      'R='+genome_fasta,
                      'O='+dict_file]
        print ' '.join(pg_command)
        process = subprocess.Popen(pg_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(err)
            sys.exit()
    else:
        print '.dict file for genome exists. OK!'
    
    # Using samtools prepare genomic fasta index to be used in GATK 
   
    fai_file = genome_fasta+'.fai'
    if not os.path.isfile(fai_file):
        pg_command = ['samtools','faidx',genome_fasta]
        print ' '.join(pg_command)
        process = subprocess.Popen(pg_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(err)
            sys.exit()

    else:
        print '.fai file for genome exists. OK!'   
    print '-----'
        
def prepare_bam(bam_file, path_to_picard):

    # Using picard tools prepare BAM file for usage in GATK: add read group and index.
    
    print '-----prepare-bam'
    assert bam_file.endswith('.bam')
    if not os.path.isfile(bam_file[:-1]+'i'):
        sample_name = bam_file.split('.')[0]
        rg_bam_file = bam_file[:-4]+'.rg.bam'
        rg_command = ['java','-Xmx1g','-jar',path_to_picard,
                      'AddOrReplaceReadGroups','RGLB=1','RGPL=illumina', 
                      'RGPU=1','RGSM='+sample_name,'CREATE_INDEX=true',
                      'I='+bam_file,
                      'O='+rg_bam_file]
        print ' '.join(rg_command)
        process = subprocess.Popen(rg_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(err)
            sys.exit()
        print 'mv %s %s' % (rg_bam_file,bam_file)  
        os.rename(rg_bam_file,bam_file) # write to initial file
        print 'mv %s %s' % (rg_bam_file[:-1]+'i',bam_file[:-1]+'i')
        os.rename(rg_bam_file[:-1]+'i',bam_file[:-1]+'i') # move index
    else:
        print 'bam file is prepared for GATK. If it is not, remove associated *bai file and repeat.'     
    print '-----'
    

def run_haplotypecaller(bam_file, path_to_gatk, genome_fasta, gatk_mem):

    # run variant calling with GATK HaplotypeCaller
    
    print '-----run-haplotypecaller'
    assert bam_file.endswith('.bam')
    vcf_file = bam_file[:-4]+'.hc.vcf'
    if not os.path.isfile(vcf_file):
        hc_command = ['java','-Xmx'+gatk_mem,'-jar',path_to_gatk,
                      '-T','HaplotypeCaller',
                      '-R',genome_fasta,
                      '-I',bam_file,
                      '-o',vcf_file]
        print ' '.join(hc_command)
        process = subprocess.Popen(hc_command)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #(out, err) = process.communicate()
        #sys.stdout.write(out)
        #sys.stderr.write(err)
        process.wait()
        if process.returncode != 0:
            sys.exit()
    else:
        print vcf_file + ' whole-genome vcf file exists. OK!'
    print '-----'

def select_region_variants(bam_file, bed_file, path_to_gatk, genome_fasta, stats=True):
    
    # select variants based on BED file with regions

    print '-----select-region-variants'
    assert bam_file.endswith('.bam')
    with open(bed_file) as infile: # bed sanity check
        for line in infile:
            if len(line) > 1: # non-empty lines
                assert line.startswith('chr'), "Improper chromosome name in bed file:\n%s" % line
                ll = line.split('\t')
                assert len(ll) == 3, "Incorrect separation or column number in bed file:\n%s" % line
                assert int(ll[1]) < int(ll[2]), "Start coordinate is larger than the end in bed file:\n%s" % line                 

    vcf_file = bam_file[:-4]+'.hc.vcf'
    reg_vcf_file = bam_file[:-4]+'.reg.hc.vcf'
    if not os.path.isfile(reg_vcf_file):
        sv_command = ['java','-Xmx1g','-jar',path_to_gatk,
                      '-T','SelectVariants',
                      '-R',genome_fasta,
                      '--variant',vcf_file,
                      '-L',bed_file,
                      '-o',reg_vcf_file]
        print ' '.join(sv_command)
        process = subprocess.Popen(sv_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(err)
            sys.exit()
    else:
        print reg_vcf_file + ' region vcf file exists. OK!'
    
    # calculate statistics of the resuling file 

    reg_vcf_stat_file = bam_file[:-4]+'.reg.hc.txt'
    if not os.path.isfile(reg_vcf_stat_file) and stats:
        ve_command = ['java','-Xmx1g','-jar',path_to_gatk,
                      '-T','VariantEval',
                      '-R',genome_fasta,
                      '--eval',reg_vcf_file,
                      '-o',reg_vcf_stat_file,
                      '-noEV', '-noST',
                      '-EV', 'CountVariants',
                      '-EV', 'TiTvVariantEvaluator',
                      '-EV', 'IndelSummary',
                      '-EV', 'MultiallelicSummary',
                      '-EV', 'ValidationReport',
                      '-ST', 'Sample',
                      '-l', 'INFO']
        print ' '.join(ve_command)
        process = subprocess.Popen(ve_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(err)
            sys.exit()
    else:
        print reg_vcf_stat_file + ' region variant stat file exists. OK!'
    print '-----'

def annotate_region_variants(bam_file, path_to_snpEff, genome_snpEff):
    
    print '-----annotate-region-variants'
    
    # annotate variants for positions residing inside target regions
    
    vcf_file = bam_file[:-3]+'reg.hc.vcf'
    ann_vcf_file = bam_file[:-3]+'ann.vcf'
    if not os.path.isfile(ann_vcf_file):
        se_command = ['java','-Xmx1g','-jar',path_to_snpEff,
                          '-stats',ann_vcf_file+'.csv', '-csvStats',
                          genome_snpEff, vcf_file]
        with open(ann_vcf_file,'w') as out_vcf:   
            print ' '.join(se_command)
            process = subprocess.Popen(se_command, stdout=out_vcf, stderr=subprocess.PIPE)
            (out, err) = process.communicate()
            if process.returncode != 0:
                sys.stderr.write(err)
                sys.exit()
    else:
        print ann_vcf_file + ' annotated vcf file exists. OK!'

    
    print '-----'
    
def calc_annot_stats(bam_file,reg_file,path_to_snpEff,genome_snpEff):

    print '-----calc-annot-stats'

    # create bed with positions inside target regions divided into 1bp chunks

    pos_file = bam_file[:-3] + 'pos.bed'
    regpos_file = bam_file[:-3] + 'regpos.bed'
    if not os.path.isfile(regpos_file):
        bi_command = ['bedtools','intersect', '-a', pos_file, '-b', reg_file]
        print ' '.join(bi_command) + ' > ' + regpos_file + ' # Splitting into 1bp chunks'
        process = subprocess.Popen(bi_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(err)
            sys.exit()
        with open(regpos_file, 'w') as outfile:
            for line in out.split('\n'):
                ll = line.split()
                if len(ll) > 2:
                    start = int(ll[1])
                    end = int(ll[2])
                    while start < end:
                        outfile.write( '%s\t%d\t%d\n' % (ll[0],start,start+1) )
                        start += 1  
    else:
        print regpos_file + ' file with positions inside target regions exists. OK!'

    # create regpos count file

    count_file = bam_file[:-3] + 'regpos.count'
    if not os.path.isfile(count_file + '.txt'):
        count_command = ['java','-Xmx1g','-jar',path_to_snpEff,'count',
                          '-n',count_file, genome_snpEff, regpos_file]
        print ' '.join(count_command)
        process = subprocess.Popen(count_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(err)
            sys.exit()
    else:
        print count_file + ' feature count files for positions inside target regions exists. OK!'
    
    # create variant density file
    
    dens_file = bam_file[:-3] + 'ann.vcf.dens.txt'
    if not os.path.isfile(dens_file):
        
        # parse variant count csv
        
        var_count = dict()
        with open(bam_file[:-3]+'ann.vcf.csv') as infile:
            i = 0
            for line in infile:
                if i == 1:
                    if line.startswith('#'):
                        break
                    ll = line.split(' , ')
                    if len(ll) == 3 and ll[0] != 'Type': # parse non-empty lines and skip header
                        var_count[ll[0]] = int(ll[1])
                elif line.startswith('# Count by genomic region'):
                    i = 1
        
        # parse base count file to obtain matching gene features
        # this step is based on the assupmtion that overlapping variant types coinside with overlapping gene annotations

        bp_count = {'SPLICE_SITE_REGION' : 0}
        with open(count_file+'.summary.txt') as infile:
            for line in infile:
                ll = line.split()
                if ll[0] in ('Downstream','Exon','Intergenic','Intron','Upstream'):
                    bp_count[ll[0].upper()] = int(ll[1])
                elif ll[0] == 'Utr3prime':
                    bp_count['UTR_3_PRIME'] = int(ll[1])
                elif ll[0] == 'Utr5prime':
                    bp_count['UTR_5_PRIME'] = int(ll[1])
                elif ll[0].startswith('SpliceSite'): # Donor, Acceptor and Region
                    bp_count['SPLICE_SITE_REGION'] += int(ll[1])

        # calculate density
        
        with open(dens_file, 'w') as out:
            print dens_file + ' variant density file writing.'
            features = bp_count.keys()
            features.sort()
            for f in features:
                if f in var_count.keys():
                    out.write( '%s\t%d\n' % (f,bp_count[f]/var_count[f]) ) # int rounding ok?
                else: # no variants in feature
                    out.write( '%s\tNA\n' % (f) )
             
        # clean up
        #print 'rm ' + regpos_file
        #os.remove(regpos_file)
    else:
        print dens_file + ' file with variant densities exists. OK!'

    print '-----'

def main(config_file):
    
    parser = ConfigParser.ConfigParser()
    parser.read(config_file)
    prepare_genome(parser.get('VC','genome_fa'), parser.get('VC','path_to_picard'))
    prepare_bam(parser.get('VC','bam_file'), parser.get('VC','path_to_picard'))
    run_haplotypecaller(parser.get('VC','bam_file'), parser.get('VC','path_to_gatk'),
                        parser.get('VC','genome_fa'), parser.get('VC','gatk_mem'))
    select_region_variants(parser.get('VC','bam_file'), parser.get('VC','reg_bed'),
                           parser.get('VC','path_to_gatk'), parser.get('VC','genome_fa'), stats=True) 
        
    if parser.get('VA','path_to_snpEff'):
        annotate_region_variants(parser.get('VC','bam_file'),
                                 parser.get('VA','path_to_snpEff'), parser.get('VA','genome_snpEff'))
        calc_annot_stats(parser.get('VC','bam_file'), parser.get('VC','reg_bed'),
                         parser.get('VA','path_to_snpEff'), parser.get('VA','genome_snpEff'))

if __name__ == '__main__':
    main(sys.argv[1])
