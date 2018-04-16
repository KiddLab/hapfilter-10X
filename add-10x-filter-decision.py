import sys
import gzip
import hapfilt

from optparse import OptionParser

###############################################################################
USAGE = """
add-10x-filter-decision.py --in <10X SNV info talbe>

Assumes/only deails with biallelic SNVs
Has some built in weak filters for quality.
will write output to infile.filter 



"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inFileName', help = '10X input file name')

(options, args) = parser.parse_args()


if options.inFileName is None:
    parser.error('in file not given')
    
###############################################################################
def parse_header(line):
    d = {}
    if line[0] == '#':
        line = line[1:]
    line = line.split()
    for i in range(len(line)):
        d[line[i]] = i
    return d



###############################################################################

of = options.inFileName + '.filter'
inFile = open(options.inFileName,'r')

print options.inFileName
print of

outFile = open(of,'w')
for line in inFile:
    line = line.rstrip()
    if line[0] == '#':
#        print 'header!'
        headerDict = parse_header(line)
#        print headerDict
        headerLine = line + '\t10xdecision\n'
        outFile.write(headerLine)
    else:
        line = line.rstrip()
        line = line.split()

        counts = [[0,0,0,0.0],[0,0,0,0.0],[0,0,0,0.0],[0.0,0.0]]
        counts[0][0] = int(line[headerDict['hap1Ref']])
        counts[0][1] = int(line[headerDict['hap1Alt']])
        counts[0][2] = int(line[headerDict['hap1Other']])        
        counts[0][3] = float(line[headerDict['hap1Freq']])        

        counts[1][0] = int(line[headerDict['hap2Ref']])
        counts[1][1] = int(line[headerDict['hap2Alt']])
        counts[1][2] = int(line[headerDict['hap2Other']])        
        counts[1][3] = float(line[headerDict['hap2Freq']])        

        counts[2][0] = int(line[headerDict['hapUnkRef']])
        counts[2][1] = int(line[headerDict['hapUnkAlt']])
        counts[2][2] = int(line[headerDict['hapUnkOther']])        
        counts[2][3] = float(line[headerDict['hapUnkFreq']])        

        counts[3][0] = float(line[headerDict['minFreq']])        
        counts[3][1] = float(line[headerDict['maxFreq']])        


        
        decision = hapfilt.get_decision(counts)
        line.append(decision)
        line = '\t'.join(line) + '\n'
        outFile.write(line)

inFile.close()
outFile.close()

