#hapfilt.py
import sys
import pysam
#######################################################################
def count_site(bamFileName,chrom,pos,ref,alt):
   # counts will be a list of list for
   # hap1, hap2, haplotype not assigned
   # and then for each haplotype
   # ref counts, alt counts, other allele, alt/(ref+alt)
   # then last two are min and max of fraction for hap1 and hap2
    counts = [[0,0,0,0.0],[0,0,0,0.0],[0,0,0,0.0],[0.0,0.0]]
    samfile = pysam.AlignmentFile(bamFileName, 'rb')
    # get  reads overlapping segment
    for read in samfile.fetch(chrom,pos-1,pos):
        # some simple filters, not being very strict at all...
        if read.is_duplicate is True:
            continue
        if read.is_qcfail is True:
            continue    
        if read.mapping_quality <= 10:
            continue    
        
        aligned_pairs = read.get_aligned_pairs(with_seq=False) 
        queryAlignedBase = '?' # incase we don't have an aligned base, probably not needed
        numAligned = 0
        for i in aligned_pairs:
            if i[1] == pos-1:
                if i[0] == None: # deletion
                    continue
                queryAlignedBase = read.query_sequence[i[0]]
                numAligned += 1
        if numAligned != 1:  #might see multiple aligned in case of indel
            continue
        if queryAlignedBase == '?' or queryAlignedBase == 'N':  # no alignment...
            continue
        if queryAlignedBase == ref:
            allele_i = 0
        elif queryAlignedBase == alt:
            allele_i = 1
        else:
            allele_i = 2
        if read.has_tag('HP'):
            hapTag = read.get_tag('HP')
            hap_num = hapTag - 1
        else:
            hap_num = 2 # no hap
        counts[hap_num][allele_i] += 1
    samfile.close()
    
    for i in range(3): # get the alt read fraction
        t = counts[i][0]+counts[i][1]
        if t == 0:
            counts[i][3] = 0
        else:
            counts[i][3] = float(counts[i][1])/t    
    

    counts[3][0] = min(counts[0][3],counts[1][3])
    counts[3][1] = max(counts[0][3],counts[1][3])
    
    return counts
#######################################################################
def get_decision(counts):
    decision = '?'
    minFreq = counts[3][0]
    maxFreq = counts[3][1]
    
    if maxFreq >= 0.90:
       decision = 'fail10Xhet'
       return decision
    if maxFreq >= 0.50:
       decision = 'fail10Xmajorhap'
       return decision
    if minFreq >= 0.01:
       decision = 'fail10Xbothhaps'
       return decision
    if maxFreq == 0.0:
        decision = 'PASS10Xnoinfo'
        return decision
    decision = 'PASS10Xok'

    
    return decision    




