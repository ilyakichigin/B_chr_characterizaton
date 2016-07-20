#!/usr/bin/env python

# From bam file create annotated VCF callset for target genomic regions and statistics. See .conf for details.
# Requires: GATK (tested on v.3.3.0), picard-tools (tested on v.1.125)
# Requires in $PATH: samtools (in $PATH, tested on v.0.1.19), bedtools
# Usage: vca_reg.py vca_reg.conf

import subprocess
import sys
import ConfigParser
import os

def prepare_genome(genome_fasta, path_to_picard):

    # Using picard tools prepare genomic fasta to be used in GATK    
    
    print '-----prepare-genome'
    assert os.path.isfile(genome_fasta)
    assert os.path.isfile(path_to_picard)
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
    assert os.path.isfile(bam_file)
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
        print '.bam file is prepared for GATK. If it is not, remove associated *bai file and repeat.'     
    print '-----'
    

def run_haplotypecaller(bam_file, path_to_gatk, genome_fasta, gatk_mem):

    # run variant calling with GATK HaplotypeCaller
    
    print '-----run-haplotypecaller'
    assert os.path.isfile(path_to_gatk)
    vcf_file = bam_file[:-3]+'hc.vcf'
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

def select_region_variants(bam_file, out_prefix, bed_file, path_to_gatk, genome_fasta, stats=True, verbose=True):
    
    # select variants based on BED file with regions

    if verbose:
        print '-----select-region-variants'

    # bed sanity check
    assert os.path.isfile(bed_file)
    with open(bed_file) as infile: 
        for line in infile:
            if len(line) > 1: # non-empty lines
                # this assert does not work for scaffold assemblies
                #assert line.startswith('chr'), "Improper chromosome name in bed file:\n%s" % line
                ll = line.split('\t') # otherwise, bedtools won't recognize the file
                assert len(ll) == 3, "Incorrect separation or column number in bed file:\n%s" % line
                assert int(ll[1]) < int(ll[2]), "Start coordinate is larger than the end in bed file:\n%s" % line

    vcf_file = bam_file[:-3]+'hc.vcf'
    reg_vcf_file = out_prefix+'hc.reg.vcf'
    reghz_vcf_file = out_prefix+'hc.reghz.vcf'
    reg_vcf_stat_file = out_prefix+'hc.reg.txt'
    reghz_vcf_stat_file = out_prefix+'hc.reghz.txt'

    # select region variants
    if not os.path.isfile(reg_vcf_file):
        sv_command = ['java','-Xmx1g','-jar',path_to_gatk,
                      '-T','SelectVariants',
                      '-R',genome_fasta,
                      '--variant',vcf_file,
                      '-L',bed_file,
                      '-o',reg_vcf_file]
        if verbose:
            print ' '.join(sv_command)
        process = subprocess.Popen(sv_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(err)
            sys.exit()
    else:
        print reg_vcf_file + ' region vcf file exists. OK!'

    # select region heterozygous variants
    sample_name = bam_file.split('.')[0]
    if not os.path.isfile(reghz_vcf_file):
        sv_command = ['java','-Xmx1g','-jar',path_to_gatk,
                      '-T','SelectVariants',
                      '-R',genome_fasta,
                      '--variant',vcf_file,
                      '-select','vc.getGenotype(\'%s\').isHet()' % (sample_name),
                      '-L',bed_file,
                      '-o',reghz_vcf_file]
        if verbose:
            print ' '.join(sv_command)
        process = subprocess.Popen(sv_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(err)
            sys.exit()
    else:
        print reghz_vcf_file + ' region heterozygous vcf file exists. OK!'

    # calculate statistics of the resuling region vcfs 
    for (vcf_file, stat_file) in [(reg_vcf_file,reg_vcf_stat_file),
                                  (reghz_vcf_file,reghz_vcf_stat_file)]:
        if not os.path.isfile(stat_file) and stats:
            ve_command = ['java', '-Xmx1g','-jar',path_to_gatk,
                          '-T', 'VariantEval',
                          '-R', genome_fasta,
                          '--eval', vcf_file,
                          '-o', stat_file,
                          '-noEV', '-noST',
                          '-EV', 'CountVariants',
                          '-EV', 'TiTvVariantEvaluator',
                          '-EV', 'IndelSummary',
                          '-EV', 'MultiallelicSummary',
                          '-EV', 'ValidationReport',
                          '-ST', 'Sample',
                          '-l', 'INFO']
            if verbose:
                print ' '.join(ve_command)
            process = subprocess.Popen(ve_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out, err) = process.communicate()
            if process.returncode != 0:
                sys.stderr.write(err)
                sys.exit()
        else:
            print stat_file + ' region variant stat file exists. OK!'
    if verbose:
        print '-----'

def generate_regpos_files(bam_file, reg_file, genome_fasta, verbose = True):

    if verbose:
        print '-----generate-regpos-files'
    assert os.path.isfile(reg_file)
    # input
    pos_file = bam_file[:-3] + 'pos.bed'
    # outputs
    regpos_file = reg_file[:-3] + 'regpos.bed'
    regpos_split_file = reg_file[:-3] + 'regpos.split.bed'

    # create bed with positions inside target and divide it into 1bp chunks
    if not os.path.isfile(regpos_file) and not os.path.isfile(regpos_split_file):

        # get chromosome order from reference
        chr_order = []
        with open(genome_fasta + '.fai') as infile:
            for line in infile:
                chr_name = line.split()[0]
                chr_order.append(chr_name)
        
        # intersect regs and pos
        bi_command = ['bedtools','intersect', '-a', pos_file, '-b', reg_file]
        if verbose:
            print '%s > %s \n%s - splitted into 1bp chunks' % (' '.join(bi_command), regpos_file, regpos_split_file )
        process = subprocess.Popen(bi_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(err)
            sys.exit()

        chr_pos = dict()
        for line in out.split('\n'):
            ll = line.split()
            if len(ll) > 2:
                if ll[0] not in chr_pos.keys():
                    chr_pos[ll[0]] = [ll]
                else:
                    chr_pos[ll[0]].append(ll)
        
        with open(regpos_file, 'w') as of, open(regpos_split_file, 'w') as splof:
            for chr_name in chr_order:
                if chr_name in chr_pos:
                    for pos in chr_pos[chr_name]:
                        of.write('\t'.join(pos) + '\n')
                        start = int(pos[1])
                        end = int(pos[2])
                        while start < end:
                            splof.write( '%s\t%d\t%d\n' % (pos[0],start,start+1) )
                            start += 1  
    else:
        print regpos_file + ' file with positions inside target regions and its splitted version exist. OK!'
    if verbose:
        print '-----'

def generate_alt_fasta(bam_file, reg_file, path_to_gatk, genome_fasta): 

    print '-----generate-alt-fasta'
    # inputs
    reg_vcf_file = bam_file[:-3]+'hc.reg.vcf'
    regpos_bed = reg_file[:-3] + 'regpos.bed'
    sample_name = bam_file.split('.')[0]
    # output
    alt_fasta_file = bam_file[:-3]+'alt.fa'
    
    if not os.path.isfile(alt_fasta_file):

        # generate alternate fasta with GATK
        af_command = ['java', '-Xmx1g', '-jar', path_to_gatk,
                      '-T', 'FastaAlternateReferenceMaker',
                      '-R', genome_fasta,
                      '-L', regpos_bed,
                      '-V', reg_vcf_file, 
                      '--use_IUPAC_sample', sample_name]
        print ' '.join(af_command)
        process = subprocess.Popen(af_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(err)
            sys.exit()

        # add reference coordinates to alternate fasta
        regpos = []
        with open(regpos_bed) as infile:
            for line in infile:
                regpos.append('_'.join(line.split()))
        with open(alt_fasta_file, 'w') as outfile:
            for line in out.split('\n'):
                if line.startswith('>'):
                    outfile.write('>'+regpos.pop(0))
                else:
                    outfile.write(line)
                outfile.write('\n')

    else:
        print alt_fasta_file + ' alternate fasta file for positions covered with reads extits. OK!'
    print '-----'

def annotate_region_variants(prefix, path_to_snpEff, genome_snpEff, verbose=True):
    
    if verbose:
        print '-----annotate-region-variants'
    assert os.path.isfile(path_to_snpEff) 
    # input
    vcf_file = prefix+'hc.reg.vcf'
    # output
    ann_vcf_file = prefix+'hc.reg.ann.vcf'

    # annotate variants for positions residing inside target regions
    if not os.path.isfile(ann_vcf_file):
        se_command = ['java','-Xmx4g','-jar',path_to_snpEff,
                          '-stats',ann_vcf_file+'.csv', '-csvStats',
                          genome_snpEff, vcf_file]
        with open(ann_vcf_file,'w') as out_vcf:   
            if verbose:
                print ' '.join(se_command) + ' > ' + ann_vcf_file
            process = subprocess.Popen(se_command, stdout=out_vcf, stderr=subprocess.PIPE)
            (out, err) = process.communicate()
            if process.returncode != 0:
                sys.stderr.write(err)
                sys.exit()
    else:
        print ann_vcf_file + ' annotated region vcf file exists. OK!'

    # input
    vcf_file = prefix+'hc.reghz.vcf'
    # output
    ann_vcf_file = prefix+'hc.reghz.ann.vcf'

    # annotate heterozygous variants for positions residing inside target regions
    if not os.path.isfile(ann_vcf_file):
        se_command = ['java','-Xmx4g','-jar',path_to_snpEff,
                          '-stats',ann_vcf_file+'.csv', '-csvStats',
                          genome_snpEff, vcf_file]
        with open(ann_vcf_file,'w') as out_vcf:   
            if verbose:
                print ' '.join(se_command) + ' > ' + ann_vcf_file
            process = subprocess.Popen(se_command, stdout=out_vcf, stderr=subprocess.PIPE)
            (out, err) = process.communicate()
            if process.returncode != 0:
                sys.stderr.write(err)
                sys.exit()
    else:
        print ann_vcf_file + ' annotated heterozygous region vcf file exists. OK!'
    
    if verbose:
        print '-----'
    
def calc_annot_stats(vcf_prefix, reg_file, path_to_snpEff, genome_snpEff, to_file = True):

    # calculate statistics of the annotated file

    if to_file:
        print '-----calc-annot-stats'

    # inputs
    regpos_split_file = reg_file[:-3] + 'regpos.split.bed'
    reg_ann_vcf_csv_file = vcf_prefix + 'hc.reg.ann.vcf.csv'
    reghz_ann_vcf_csv_file = vcf_prefix + 'hc.reghz.ann.vcf.csv'
    # outputs
    count_base = reg_file[:-3] + 'regpos.count'
    count_file = count_base + '.txt'
    count_summary_file = count_base + '.summary.txt'
    if to_file:
        dens_file = vcf_prefix + 'hc.reg.ann.dens.txt'
    else:
        dens_file = 'NA' # nonsence name for check pass


    # create regpos count file
    if not os.path.isfile(count_file):
        if os.path.isfile(dens_file):
            os.remove(dens_file)
        count_command = ['java','-Xmx4g','-jar',path_to_snpEff,'count',
                          '-n',count_base, genome_snpEff, regpos_split_file]
        if to_file:
            print ' '.join(count_command)
        process = subprocess.Popen(count_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out, err) = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(err)
            sys.exit()
    else:
        print count_file + ' feature count files for positions inside target regions exists. OK!'
    
    # create variant density file
    if not os.path.isfile(dens_file):       
        # parse variant count csv's
        # region variants
        var_count = {'TOTAL':0,'DOWNSTREAM':0,'EXON':0,'INTERGENIC':0,'INTRON':0,
                    'SPLICE_SITE_REGION':0,'UPSTREAM':0,'UTR_3_PRIME':0,'UTR_5_PRIME':0}
        with open(reg_ann_vcf_csv_file) as infile:
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
                elif line.startswith('Number_of_variants_processed'):
                    var_count['TOTAL'] = int(line.split(' , ')[1])
        # region het variants
        varhz_count = {'TOTAL':0,'DOWNSTREAM':0,'EXON':0,'INTERGENIC':0,'INTRON':0,
                    'SPLICE_SITE_REGION':0,'UPSTREAM':0,'UTR_3_PRIME':0,'UTR_5_PRIME':0}
        with open(reghz_ann_vcf_csv_file) as infile:
            i = 0
            for line in infile:
                if i == 1:
                    if line.startswith('#'):
                        break
                    ll = line.split(' , ')
                    if len(ll) == 3 and ll[0] != 'Type': # parse non-empty lines and skip header
                        varhz_count[ll[0]] = int(ll[1])
                elif line.startswith('# Count by genomic region'):
                    i = 1
                elif line.startswith('Number_of_variants_processed'):
                    varhz_count['TOTAL'] = int(line.split(' , ')[1])
        
        # parse base count file to obtain matching gene features
        # this step is based on the assupmtion that overlapping variant types coinside with overlapping gene annotations
        bp_count = {'DOWNSTREAM':0,'EXON':0,'INTERGENIC':0,'INTRON':0,
                            'SPLICE_SITE_REGION':0,'UPSTREAM':0,'UTR_3_PRIME':0,'UTR_5_PRIME':0}
        with open(count_summary_file) as infile:
            for line in infile:
                ll = line.split()
                if ll[0] in ('Downstream','Exon','Intergenic','Intron','Upstream'):
                    bp_count[ll[0].upper()] = int(ll[1])
                elif ll[0] == 'Utr3prime':
                    bp_count['UTR_3_PRIME'] = int(ll[1])
                elif ll[0] == 'Utr5prime':
                    bp_count['UTR_5_PRIME'] = int(ll[1])
                elif ll[0] == 'Chromosome':
                    bp_count['TOTAL'] = int(ll[1])
                elif ll[0].startswith('SpliceSite'): # Donor, Acceptor and Region
                    bp_count['SPLICE_SITE_REGION'] += int(ll[1])

        # calculate densities
        if to_file:
            with open(dens_file, 'w') as out:
                out.write( '#Feature\tbp_covered\tvariants\tbp/var\tvariants_hz\tbp/var_hz\n' )
                features = ['TOTAL','DOWNSTREAM','EXON','INTERGENIC','INTRON',
                            'SPLICE_SITE_REGION','UPSTREAM','UTR_3_PRIME','UTR_5_PRIME']
                for f in features:
                    if var_count[f] > 0: # variants in region 
                        if varhz_count[f] > 0: # heterozygous variants in region 
                            out.write( '%s\t%d\t%d\t%d\t%d\t%d\n' % (f,bp_count[f],var_count[f],bp_count[f]/var_count[f],varhz_count[f],bp_count[f]/varhz_count[f]) ) # note int rounding for density calculation
                        else: # no heterozygous variants in feature
                            out.write( '%s\t%d\t%d\t%d\t0\tNA\n' % (f,bp_count[f],var_count[f],bp_count[f]/var_count[f]) ) # note int rounding for density calculation
                    elif f in bp_count.keys(): # no variants in feature
                        out.write( '%s\t%d\t0\tNA\t0\tNA\n' % (f,bp_count[f]) )
                    else: # no feature in reads for positions
                        out.write( '%s\t0\t0\tNA\t0\tNA\n' % (f) )
        else:
            out = []
            features = ['TOTAL','INTERGENIC','DOWNSTREAM','UPSTREAM','INTRON','EXON',
                        'SPLICE_SITE_REGION','UTR_3_PRIME','UTR_5_PRIME']
            for f in features:
                if f in bp_count.keys():
                    out.append((f + '_bp',bp_count[f]))
                else:
                    out.append((f + '_bp',0))
                if var_count[f] > 0: # variants in region 
                    out.append((f + '_var',var_count[f]))
                    out.append((f+'_dens',bp_count[f]/var_count[f]))
                else:
                    out.append((f + '_var',0))
                    out.append((f+'_dens','NA'))
                if varhz_count[f] > 0:
                    out.append((f + '_varhz',varhz_count[f]))
                    out.append((f+'_denshz',bp_count[f]/varhz_count[f]))
                else:
                    out.append((f + '_varhz',0))
                    out.append((f+'_denshz','NA'))

            return out

    else:
        print dens_file + ' file with variant densities exists. OK!'
    if(to_file):
        print '-----'

def genes_in_reg(bam_file, reg_file):
    
    # based on previously created files, add gene names and Ensmebl IDs to region BED file.
    
    print '-----genes_in_reg'
    #inputs
    regpos_count_file = reg_file[:-3] + 'regpos.count.txt'
    vcf_genes_file = bam_file[:-3] + 'hc.reg.ann.vcf.genes.txt'
    #output
    genereg_file = reg_file[:-3] + 'genes.bed'

    if not os.path.isfile(genereg_file):
        # regions to dict
        regs = []
        with open(reg_file) as f:
            for line in f:
                ll = line.split()
                if len(ll) > 0:
                    assert len(ll) == 3 # Region file is clean BED
                    line = '\t'.join(ll) + '\n' # tab-delimited BED from space-delimited BED
                    regs.append(line)

        # extract Ensembl gene IDs from regpos count file 
        gene_ids = []
        with open(regpos_count_file) as f:
            next(f) # skip header
            for line in f:
                ll = line.split()
                if len(ll) == 6:
                    if ll[3].split(';')[0] == 'Gene':
                        gene_id = ll[3].split(';')[1]
                        chrom = ll[0]
                        if chrom[0].isdigit() or chrom == 'X' or chrom == 'Y': # assume that in reg.bed chromosomes have 'chr' prefix
                            chrom = 'chr' + chrom
                        start = ll[1]
                        end = ll[2]
                        gene_ids.append('\t'.join([chrom,start,end,gene_id]))
                        # match gene coordinates with region
        
        # match gene ids with individual regions based on genomic coordinates.
        gene_bed_file = 'tmp.ensGene.bed'
        with open(gene_bed_file, 'w') as f:
            f.write('\n'.join(gene_ids))
        regs_ensGene = dict()
        reg_bed_file = 'tmp.reg.bed'
        for reg in regs:
            with open(reg_bed_file, 'w') as f:
                f.write(reg)
            bi_command = ['bedtools','intersect', '-a', gene_bed_file, '-b', reg_bed_file]
            #print ' '.join(bi_command) 
            process = subprocess.Popen(bi_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out, err) = process.communicate()
            if process.returncode != 0:
                sys.stderr.write(err)
                sys.exit()
            reg_eg = [line.split()[3] for line in out.split('\n')[:-1]]
            regs_ensGene[reg] = reg_eg
            os.remove(reg_bed_file)
        os.remove(gene_bed_file)

        # match conventional gene names with Ensembl gene IDs
        gene_names = dict()
        with open(vcf_genes_file) as f:
            for line in f:
                ll = line.split()
                gene_names[ll[1]] = ll[0]

        # get gene names from Ensembl IDs for individual regions and write to file
        with open(genereg_file, 'w') as out:
            
            print 'Adding genes to region BED %s. Writing to %s.' % (reg_file,genereg_file)
            for reg in regs:
                ensGenes = regs_ensGene[reg]
                if len(ensGenes) > 0:
                    namedGenes = []
                    for ensGene in ensGenes:
                        if ensGene in gene_names.keys():
                           namedGenes.append(gene_names[ensGene])
                        else:
                            namedGenes.append(ensGene)
                else: 
                    ensGenes = ['NA']
                    namedGenes = ['NA']
                out.write( '\t'.join([reg.strip('\n'), ';'.join(ensGenes), ';'.join(namedGenes)]) + '\n' )
    
    else:
         print genereg_file + ' region BED file with genes added exists. OK!'

    print '-----'

def tmp_cleanup(prefix):
    filelist = [ file for file in os.listdir(".") if file.startswith(prefix) ]
    for file in filelist:
        os.remove(file)

def reg_calc_annot_stats(bam_file, reg_file, path_to_gatk, path_to_snpEff, genome_snpEff, genome_fasta):

    # calculate annotation statistics (as in calc_annot_stats()) independently for each region
    # add statistics to bed with genes (output of genes_in_reg())

    print '----reg_calc_annot_stats'
    # input
    genes_file = reg_file[:-3] + 'genes.bed'
    # output
    annot_file = reg_file[:-3] + 'ann.genes.txt'

    if not os.path.isfile(annot_file):
        with open(genes_file) as f, open(annot_file, 'w') as o:

            tmp_prefix = 'tmp.one.'
            i = 0
            header = True
            o.write('chr\tstart\tend\tensGene\tgene')
            for line in f: # same output for each line
                tmp_cleanup(tmp_prefix)
                i += 1
                sys.stdout.write('Processing region %s ' % (i))
                # write one-region BED
                one_reg_file = tmp_prefix + reg_file
                with open(one_reg_file, 'w') as of:
                    bed = '\t'.join(line.split()[:3])
                    of.write(bed + '\n')
                sys.stdout.write('Selecting variants ')
                sys.stdout.flush()
                select_region_variants(bam_file, tmp_prefix, one_reg_file, path_to_gatk, genome_fasta, stats=True, verbose=False)
                sys.stdout.write('Generating regpos ')
                sys.stdout.flush()
                generate_regpos_files(bam_file, one_reg_file, genome_fasta, verbose = False)
                sys.stdout.write('Annotating variants ')
                sys.stdout.flush()
                annotate_region_variants(tmp_prefix, path_to_snpEff, genome_snpEff, verbose=False)
                sys.stdout.write('Calculating stats ')
                sys.stdout.flush()
                stats = calc_annot_stats(tmp_prefix, one_reg_file, path_to_snpEff, genome_snpEff, to_file = False)
                #print stats
                if header:
                    for j in range(len(stats)):
                        o.write('\t'+stats[j][0])
                        header = False
                    o.write('\n')
                o.write(line[:-1])
                for j in range(len(stats)):
                    o.write('\t'+str(stats[j][1]))
                o.write('\n')
                sys.stdout.write('\n')

        tmp_cleanup(tmp_prefix)
    else:
        print annot_file + ' region BED file with genes and per region statistics on variation added exists. OK!'

    print '----'

def main(config_file):
    
    parser = ConfigParser.ConfigParser()
    parser.read(config_file)
    prepare_genome(parser.get('VC','genome_fasta'), parser.get('VC','path_to_picard'))
    prepare_bam(parser.get('VC','bam_file'), parser.get('VC','path_to_picard'))
    run_haplotypecaller(parser.get('VC','bam_file'), parser.get('VC','path_to_gatk'),
                        parser.get('VC','genome_fasta'), parser.get('VC','gatk_mem'))
    select_region_variants(parser.get('VC','bam_file'), parser.get('VC','bam_file')[:-3], 
                           parser.get('VC','reg_bed'), parser.get('VC','path_to_gatk'), parser.get('VC','genome_fasta'))
    generate_regpos_files(parser.get('VC','bam_file'), parser.get('VC','reg_bed'),
                          parser.get('VC','genome_fasta'))
    generate_alt_fasta(parser.get('VC','bam_file'), parser.get('VC','reg_bed'),
                        parser.get('VC','path_to_gatk'), parser.get('VC','genome_fasta')) 
        
    # annotate only if path to snpEff is given in conf
    if parser.get('VA','path_to_snpEff'):
        annotate_region_variants(parser.get('VC','bam_file')[:-3],
                                 parser.get('VA','path_to_snpEff'), parser.get('VA','genome_snpEff'))
        calc_annot_stats(parser.get('VC','bam_file')[:-3], 
                         parser.get('VC','reg_bed'), parser.get('VA','path_to_snpEff'), 
                         parser.get('VA','genome_snpEff'))
        genes_in_reg(parser.get('VC','bam_file'), parser.get('VC','reg_bed'))
        reg_calc_annot_stats(parser.get('VC','bam_file'), parser.get('VC','reg_bed'),
                         parser.get('VC','path_to_gatk'), parser.get('VA','path_to_snpEff'), 
                         parser.get('VA','genome_snpEff'), parser.get('VC','genome_fasta'))
if __name__ == '__main__':
    main(sys.argv[1])
