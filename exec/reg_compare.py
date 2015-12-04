#!/usr/bin/env python

# Compares pairwise distances for two samples. Outputs significantly differing regions.
# Additional filtering of input regions by size and number of positions and filtering of
# comparison results by pd ratio can be adjusted in function parameters.
# Usage: reg_compare.py regs1.tsv regs2.tsv

import sys

def read_regions(filename):

    regs = []
    with open(filename) as rf:
        rf.next()
        for line in rf:
            regs.append(line.split('\t'))
    return regs

def find_corr_regions(reg1, regs2, min_size = 10000, min_pos = 5, min_overlap = 5000):

    corr_regs2 = []

    reg1_chr = reg1[0]
    reg1_start = int(reg1[1])
    reg1_end = int(reg1[2])
    reg1_pos = int(reg1[3])
    assert reg1_end > reg1_start

    if (reg1_end - reg1_start >= min_size) and (reg1_pos >= min_pos): # reg1 basic condition
        for reg2 in regs2:
            reg2_chr = reg2[0]
            if reg2_chr == reg1_chr:
                reg2_start = int(reg2[1])
                reg2_end = int(reg2[2])
                reg2_pos = int(reg2[3])
                if (reg2_end - reg2_start >= min_size) and (reg2_pos >= min_pos): # reg2 basic condition
                    if (reg2_end > reg1_start + min_overlap) and (reg2_start < reg1_end - min_overlap): # reg1-reg2 overlap condition
                        corr_regs2.append(reg2)

    return corr_regs2
 
def filter_corr_regions(reg1, corr_regs2, pd_ratio = 4): # reg2:reg1 pd ratio 

    filtered_matches = []
    pd_ratio = float(pd_ratio)
    
    for reg2 in corr_regs2:
        if float(reg2[5])/float(reg1[5]) >= pd_ratio:
            filtered_matches.append('reg1\t'+'\t'.join(reg1))
            for reg2 in corr_regs2:
                filtered_matches.append('reg2\t'+'\t'.join(reg2))
            filtered_matches.append('\n')

    return filtered_matches
            

if __name__ == '__main__':

    filtered_matches = []

    regs1 = read_regions(sys.argv[1])
    regs2 = read_regions(sys.argv[2])

    for reg1 in regs1:
        corr_regs2 = find_corr_regions(reg1, regs2)
        filtered_matches = filter_corr_regions(reg1, corr_regs2)
        if len(filtered_matches) > 0:
            for m in filtered_matches:
                sys.stdout.write(m)