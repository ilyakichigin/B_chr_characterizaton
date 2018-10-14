# supress warnings for rpy2-pandas interface deprecation
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pybedtools
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

in_bam = snakemake.input[0]
genome_fai = snakemake.input.genome_fai
out_pos = snakemake.output.pos
out_reg = snakemake.output.reg
sample = snakemake.params.sample

# genome sizes
genome = dict()
with open(genome_fai) as f:
    for l in f:
        ll = l.split('\t')
        try:
            genome[ll[0]] = (0, int(ll[1]))
        except:
            pass

# positions and distances between positions
in_bam = pybedtools.BedTool(in_bam)
pos_bed = in_bam.bam_to_bed().sort().merge(c=1, o='count').saveas(out_pos)
pos = pos_bed.to_dataframe()
pos['chrom'] = pos['chrom'].astype(str)
dist = pos_bed.complement(g=genome).to_dataframe()
dist['chrom'] = dist['chrom'].astype(str)
dist['log.dist'] = np.log10(dist['end'] - dist['start'])

# distance-based segmentation of genome
dnacopy = importr('DNAcopy')
cna = robjects.r['CNA'](robjects.FloatVector(dist['log.dist']),
                        robjects.StrVector(dist['chrom']), 
                        robjects.IntVector(dist['end']), 
                        data_type="logratio", sample=sample)
cna = robjects.r['smooth.CNA'](cna)
segm = robjects.r['segment'](cna, verbose=0)

# data to python, start annotation
pandas2ri.activate()
regions = robjects.pandas2ri.ri2py(segm[1])
regions['pd.mean'] = np.power(10, regions['seg.mean'])
regions = regions.rename(columns={col: col.replace('.','_') for col in regions.columns})

# expand high-density regions by 1 position at left margin
# this step is needed as left-side distances are used in segmentation
# note that we shift margins without recalculating distances
shift_regs = []
# slow iterrows
for r in regions.iterrows():
    cr = r[1]
    # shift end back 1bp if current: high pd, next: low pd
    try:
        nr = regions.iloc[r[0] + 1]
        if cr.chrom == nr.chrom and cr.seg_mean > nr.seg_mean:
            cr['loc_end'] = cr.loc_end - 1
    except:
        pass
    # shift start back 1pos if previous: high bp, current: low pd
    try:    
        pr = regions.iloc[r[0] - 1]
        if pr.chrom == cr.chrom and pr.seg_mean > cr.seg_mean:
            cr['loc_start'] = pr.loc_end
    except:
        pass
    shift_regs.append(cr)
shift_regs = pd.DataFrame(shift_regs)
shift_regs = shift_regs.rename(columns={
        'ID':'sample',
        'loc_start':'reg_start', 
        'loc_end':'reg_end',
        'num_mark':'reg_pos',
        'seg_mean':'lg_pd_mean'})

# main annotation
def segment_positions(segm):
    """
    Using positions within a segment, correct start and end of segment to start
    of first and end of last position in segment,
    recalculate number of positions (initial clustering was done on complements 
    - one more position),
    get mean coverage and position size for positions within a segment
    """
    segm_pos = pos[(pos['chrom'] == segm['chrom']) & 
                    (pos['start'] >= segm['reg_start']) & 
                    (pos['start'] <= segm['reg_end'])]
    try:
        segm['reg_start'] = segm_pos.iloc[0]['start'] # start of first pos in segment
        segm['reg_end'] = segm_pos.iloc[-1]['end'] # end of last pos in segment
        segm['reg_pos'] = segm_pos.shape[0] # nrows
        segm['pos_cov_mean'] = segm_pos['name'].mean() # coverage as 'name' column in pos df
        segm['reg_reads'] = segm_pos['name'].sum()
        segm['pos_len_mean'] = (segm_pos['end']-segm_pos['start']).mean()
        segm['pos_len_sum'] = (segm_pos['end']-segm_pos['start']).sum()
    except: # no positions in chromosome 
        segm['reg_start'] = 0
        segm['reg_end'] = segm['reg_end'] 
        segm['reg_pos'] = 0
        segm['pos_cov_mean'] = 0
        segm['reg_reads'] = 0
        segm['pos_len_mean'] = 0
        segm['pos_len_sum'] = 0
        # replacing seg_mean and pd_mean, as values produced by DNAcopy are not meaningful
        segm['pd_mean'] = segm['reg_end']
        segm['lg_pd_mean'] = np.log10(segm['reg_end'])
    segm['reg_len'] = segm['reg_end'] - segm['reg_start']
    segm['chrom_len'] = genome[segm['chrom']][1]
    segm['prop_reg_covered'] = segm['pos_len_sum'] / segm['reg_len']
    segm['prop_chrom_len'] = segm['reg_len'] / segm['chrom_len']
    return segm
# slowest part of script - apply with entire positions df filtering at each segment
corr_regs = shift_regs.apply(segment_positions, axis=1)

#output
corr_regs.to_csv(out_reg, sep="\t", index=False)