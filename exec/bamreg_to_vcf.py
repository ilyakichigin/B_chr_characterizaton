#!/usr/bin/env python


import subprocess
import sys
import argparse
import os


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    From bam file create VCF callset for genomic regions specified in region bed file. 
                    Requires: GATK (tested on v.3.3.0), picard-tools (tested on v.1.125), samtools (in $PATH, tested on v.0.1.19)
                    """
                    )
    parser.add_argument("bam_file", help="input BAM file")
    parser.add_argument("-t", help="input BED file with target regions")
    parser.add_argument("-g",help="path to reference genome fasta file (unpacked)")
    parser.add_argument("--path_to_gatk", help="path to GATK jar file")
    parser.add_argument("--path_to_picard", help="path to picard jar file")
    parser.add_argument("--gatk_mem", default='8g', help="Memory allocated to GATK HaplotypeCaller. For other commands, 1g is allocated")
    return parser.parse_args()

def prepare_genome(genome_fasta, path_to_picard):

    # Using picard tools prepare genomic fasta to be used in GATK    
    
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
    
    assert bam_file.endswith('.bam')
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
    print '-----'
    

def run_haplotypecaller(bam_file, path_to_gatk, genome_fasta, gatk_mem):

    # run variant calling with GATK HaplotypeCaller
    
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
        print 'whole-genome vcf file exists. OK!'
    print '-----'

def select_region_variants(bam_file, bed_file, path_to_gatk, genome_fasta, stats=True):
    
    # select variants based on BED file with regions

    assert bam_file.endswith('.bam')
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
        print 'region vcf file exists. OK!'
    print '-----'

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
        print 'region variant stat file exists. OK!'

    
def main(args):
    prepare_genome(args.g, args.path_to_picard)
    prepare_bam(args.bam_file, args.path_to_picard)
    run_haplotypecaller(args.bam_file,args.path_to_gatk,args.g,args.gatk_mem)
    select_region_variants(args.bam_file, args.t, args.path_to_gatk, args.g, stats=True) 
  
if __name__ == '__main__':
    main(parse_command_line_arguments())
    

