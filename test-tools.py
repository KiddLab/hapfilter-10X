import sys
import hapfilt

bamFileName = '/home/jmkidd/links/dcmb-brainsom/technical/working/20170601_10x_bam_longranger/brain.bam'

#SNVval1	1	228475850	.	A	G	chr1:228475791-228476023
chrom = '1'
pos = 228475850
ref = 'A'
alt = 'G'


chrom = '19'
pos = 4510979
ref = 'A'
alt = 'C'


#SNVval7	6	163235167	.	T
chrom = '6'
pos = 163235167
ref = 'T'
alt = 'C'
#counts = hapfilt.count_site(bamFileName,chrom,pos,ref,alt)
#print counts

#nl = [chrom,pos,ref,alt]
#for i in counts:
#    nl.extend(i)    
#nl = [str(j) for j in nl]    
#print nl



chrom = '16'
pos = 30795009
ref = 'A'
alt = 'T'
counts = hapfilt.count_site(bamFileName,chrom,pos,ref,alt)
print counts

nl = [chrom,pos,ref,alt]
for i in counts:
    nl.extend(i)    
nl = [str(j) for j in nl]    
print nl



#22      50502469        .       A       G 
print 'next'
chrom = '22'
pos = 50502469
ref = 'A'
alt = 'G'
counts = hapfilt.count_site(bamFileName,chrom,pos,ref,alt)
print counts

nl = [chrom,pos,ref,alt]
for i in counts:
    nl.extend(i)    
nl = [str(j) for j in nl]    
print nl
