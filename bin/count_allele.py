#!/usr/bin/env python

import os
import sys
import re
from varnorm.varcharkey import VarCharKey
from varnorm.VCFNormalization import VCFNormalization

def calc_freq(mask, gtIn, a):
    if mask != None:
        gt = [gtIn[i] for i in mask]
    else:
        gt = gtIn
    
    aacnt = sum(map(lambda t: t.count(a) ,gt))
    total_ac = sum(map(lambda t: len(t), gt))
    af = float(aacnt) / total_ac
    
    total_gt = map(lambda x: len(x), gt).count(2)
    homo_alt_gt = gt.count(a*2)
    homo_alt_gtf = float(homo_alt_gt) / total_gt
    het_gt = gt.count(a + '0') + gt.count('0' + a)
    het_gtf = float(het_gt) / total_gt

    return [total_ac, aacnt, af, total_gt, homo_alt_gt, homo_alt_gtf, het_gt, het_gtf]
    

def loadMask(file):
    fi = open(file, 'rb')
    mask = {'ALL':None}
    for line in fi:
        tokens = line.strip().split('\t')
        mask[tokens[0]] = map(lambda x: int(x), tokens[1:])
    
    return mask

popMap = {'ALL': 2, 'AFR': 1, 'AMR': 3, 'EAS': 4, 'EUR': 5, 'SAS': 6}

#========================================
#
#  main function
#
#=======================================

def main():

    genome = "/sc/orga/projects/PBG/KBase/download/hgdownload.cse.ucsc.edu/goldenPath/hg19/hg19.fa"
    mask = loadMask("/sc/orga/projects/PBG/KBase/download/2500Genomes_latest/parallel_processing/2kgn_p3_v5.mask")
    vcfn = VCFNormalization(None, genome)
    fi = open(sys.argv[1], 'rb')
    fo = open(sys.argv[2], 'wb')
    for l in fi:
        if l.startswith('#'):
            continue
        tokens = l.strip().split('\t')
        
        chrom = tokens[0].replace('chr', '')
        pos = tokens[1]

        vid = tokens[2] if tokens[2] != '.' else '\N'
        ref = tokens[3]
        qual = tokens[5]
        filt = tokens[6]

        start = int(pos)
        end = start + len(ref) - 1
        
        alts = tokens[4].split(',')
        fmt = tokens[8].split(':')
        gt_idx = fmt.index("GT")

        for a in range(len(alts)):
            
            gt = map(lambda x: x.split(':')[gt_idx] , tokens[9:])
            gt = map(lambda x: re.sub('/|\||\.', '', x), gt)

            alt = alts[a]
            
            start_norm = start
            end_norm = end
            ref_norm = ref
            alt_norm = alt
            
            if len(ref) > 1 or len(alt) > 1 or ref == '.' or alt == '.' or ref == '-' or alt == '-':
                chrom, start_norm, end_norm, ref_norm, alt_norm = vcfn.normAVar(chrom, start, end, ref, alt)
            
            vkey = VarCharKey.v2k(chrom, start_norm, end_norm, alt_norm)
            for m in mask:
                freqs = calc_freq(mask[m], gt, str(a+1))
                #outline = [chrom, start, end, ref, alt, '15', '1KG_p3_v5_'+m, m, '34' ]
                outline = [vkey, 15, popMap[m], m, 34 ]
                print >> fo, '\t'.join([str(x) for x in outline  + freqs])

if __name__ == "__main__":
    sys.exit(main())


