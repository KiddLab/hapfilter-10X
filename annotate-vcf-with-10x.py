import sys
import gzip
import hapfilt

from optparse import OptionParser

###############################################################################
USAGE = """
annotate-vcf-with-10x.py --vcf <vcf file>
                         --out <output file of info>
                         --bam <10X bam file to uses> 
                         --gzip <input is gzipped, output will be too>
                         --pass <only do sites marked as PASS in vcf>


Assumes/only deails with biallelic SNVs
Has some built in weak filters for quality.


"""
parser = OptionParser(USAGE)
parser.add_option('--vcf',dest='vcfFileName', help = 'name of VCF file')
parser.add_option('--out',dest='outFileName', help = 'name of output table file')
parser.add_option('--bam',dest='bamFileName', help = 'name of 10X BAM file')
parser.add_option('--gzip',dest='isGzip',  action='store_true', default = False, help = 'input/output gzipped')
parser.add_option('--pass',dest='passOnly',  action='store_true', default = False, help = 'only do pass sites')


(options, args) = parser.parse_args()


if options.vcfFileName is None:
    parser.error('vcfFileName not given')
if options.outFileName is None:
    parser.error('outFileName not given')
if options.bamFileName is None:
    parser.error('bamFileName not given')
    
###############################################################################

if options.isGzip is True:
    inFile = gzip.open(options.vcfFileName,'r')
    outFile = gzip.open(options.outFileName,'w')
else:
    inFile = open(options.vcfFileName,'r')
    outFile = open(options.outFileName,'w')

# header
nl = ['#chrom','pos','ref','alt','passValfromVCF','hap1Ref','hap1Alt','hap1Other','hap1Freq',
      'hap2Ref','hap2Alt','hap2Other','hap2Freq',
      'hapUnkRef','hapUnkAlt','hapUnkOther','hapUnkFreq',
      'minFreq','maxFreq','10xdecision']
      
nl = '\t'.join(nl) + '\n'
outFile.write(nl)
numRecords = 0              
for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    line = line.split()
    chrom = line[0]
    pos = int(line[1])
    ref = line[3]
    alt = line[4]
    passVal = line[6]
    
    if options.passOnly is True and passVal != 'PASS':
        continue
    
    # only bialelic indels
    if len(ref) != 1:
        continue
    if len(alt) != 1:
        continue
        
    counts = hapfilt.count_site(options.bamFileName,chrom,pos,ref,alt)

    nl = [chrom,pos,ref,alt,passVal]
    for i in counts:
        nl.extend(i)    
    
    decision = hapfilt.get_decision(counts)
    nl.append(decision)

    
    nl = [str(j) for j in nl]    
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)
    numRecords += 1
    if numRecords % 100 == 0:
        print 'Did %i records...'  % numRecords
    
    
inFile.close()
outFile.close()    