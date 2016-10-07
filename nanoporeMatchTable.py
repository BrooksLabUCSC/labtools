#!/usr/bin/python

import sys, os, re, getopt, argparse, textwrap

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\

This program creates a table with adapter positions in nanopore reads, based on a psl format input file.
Optionally a second alignment file of reads vs transcripts can be added in psl or sam format, the program
then adds start and end positions of the transcript alignment in the read.
If the sam format is used, you MUST give in a file with transcript sizes in the format

TranscriptID\tsize

The program faCount outputs this format.

        '''))
group = parser.add_argument_group('required arguments')
group.add_argument('--adaptalign', type=str, required=True, help="psl format alignment of ADAPTERS vs READS")
group.add_argument('--readsizes', type=str, required=True, help="faCount output for nanopore reads")
# optional flag
parser.add_argument('--txalign', type=str,  help="psl or sam format file of READS vs TRANSCRIPTS")
parser.add_argument('--txsizes', type=str,  help="faCount output for transcripts, required if using sam format")


class readCoords(object):
    """
    keeps track of read match positions for the final output table
    """
    def __init__(self, rName, size):
        self.name = rName
        self.size = size
        self.laStart = 'NA'
        self.tStart = 'NA'
        self.tEnd = 'NA'
        self.tSize = 'NA'
        self.raStart = 'NA'
        self.tID = 'NA'
        self.tStrand = 'NA'
        self.laID = 'NA'
        self.laScore = False
        self.raID = 'NA'
        self.raScore = False
    def addLA(self, laID, laStart, laScore):
        """Add left adapter, replace if a better scoring hit is found """
        if self.laID != 'NA':
            if self.laScore > laScore:
                return
            #print >>sys.stderr, "already have {} as left adapter for {}".format(self.laID, self.name)
        self.laStart = laStart
        self.laID = laID
        self.laScore = laScore
    def addRA(self, raID, raStart, raScore):
        """Add right adapter, replace if a better scoring hit is found """
        if self.raID != 'NA':
            if self.raScore > raScore:
                return
            #print >>sys.stderr, "already have {} as right adapter for {}".format(self.raID, self.name)
        self.raStart = raStart
        self.raID = raID
        self.raScore = raScore
    def addTx(self, tName, start, end, strand):
        """Add transcript, replace if a longer hit is found"""
        skip = False
        if self.tID != 'NA':
            #print >>sys.stderr, "already have {} as tx hit for {}".format(self.tID, self.name)
            if (self.tEnd - self.tStart > end - start):
                #print >>sys.stderr, "new hit smaller than old, skipping", self.tID, tName
                skip = True
        if not skip:
            self.tStart = start
            self.tEnd = end
            self.tStrand = strand
            self.tID = tName 
            if '|' in self.tID:
                self.tID = ('|').join(tName.split('|')[:2])
            try:
                self.tSize = txSizes[self.tID]
            except KeyError:
                print >>sys.stderr, "cannot find size for", self.tID
                sys.exit()
    def printRead(self, withTx):
        """Print table to stdout """
        if withTx:
	    print "{name}\t{las}\t{ts}\t{te}\t{ras}\t{ss}\t{lid}\t{rid}\t{tid}\t{tsi}\t{str}".format(name=self.name, 
                las=self.laStart, ts=self.tStart, te=self.tEnd, ras=self.raStart, ss=self.size, 
                lid=self.laID, rid=self.raID, tid=self.tID, tsi=self.tSize, str=self.tStrand)
        else:
	    print "{name}\t{las}\t{ras}\t{ss}\t{lid}\t{rid}".format(name=self.name, 
            las=self.laStart, ras=self.raStart, ss=self.size, 
            lid=self.laID, rid=self.raID)
       
class SamHit(object):
    """
    holds one sam format blat hit
    replaces if new hit covers more of the read
    """
    def __init__(self, inline):
        fields = inline.split('\t')
        #  isMinus is 0 for + and 16 for -
        isMinus = fields[1]
        cigar = fields[5]
        seq = fields[9]
        if isMinus == '16':
            self.strand = '-'
        elif isMinus == '0':
            self.strand = '+'
        self.qName = fields[0]
        self.tName = fields[2]
        self.qSize = len(seq)
        self.qStart = 0
        self.qEnd = self.qSize
        m=re.findall("(\d+)([MIDS])", cigar)
        if (m[0][1] =='S'):
            if self.strand == '-':
                self.qEnd = self.qSize - int(m[0][0])
            else: 
                self.qStart = int(m[0][0])
        if (m[-1][1] =='S'):
            if self.strand == '-':
                self.qStart = int(m[-1][0])
            else:
                self.qEnd = self.qSize - int(m[-1][0])

class BlatHit(object):
    """
    holds one psl format blat hit
    """
    def __init__(self, inline):
        fields = inline.split('\t')
        self.match, self.mismatch, self.rep, self.Ns = map(int, fields[:4])
        self.qGapCount, self.qGapBases, self.tGapCount, self.tGapBases = map(int, fields[4:8])
        self.strand = fields[8]
        self.qName = fields[9]
        self.qSize, self.qStart, self.qEnd = map(int, fields[10:13])
        self.tName = fields[13] 
        # gencode ID
        if '|' in self.tName:
            self.tName = ('|').join(self.tName.split('|')[:2])
        #self.gName = self.tName.split('|')[1]
        self.tSize, self.tStart, self.tEnd = map(int, fields[14:17])
        self.blockCount = int(fields[17])
        self.blockSizes, self.qStarts, self.tStarts = fields[18:]
        self.qCover = self.qEnd - self.qStart	# length of blat match
        self.pslScore = self.match + self.rep - self.mismatch - self.qGapCount - self.tGapCount
        self.inline = inline

def txSizeDict(infile):
    """
    Read faCount output file and return dictionary with read lengths
    """
    sdict = dict()

    with open(infile, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('total'):
                continue
            fields = line.split('\t')
            name=fields[0]
            if '|' in name:
                name = ('|').join(name.split('|')[:2])
            sdict[name] = int(fields[1])
    return sdict

    

def sizeDict(infile):
    """
    Read faCount output file and return dictionary with read objects
    """
    sdict = dict()

    with open(infile, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('total'):
                continue
            fields = line.split('\t')
            name=fields[0]
            readObj = readCoords(name, int(fields[1]))
            sdict[name] = readObj
    return sdict


args = parser.parse_args()
# Run program

if (args.txalign and not args.txsizes) or (args.txsizes and not args.txalign):
    print >>sys.stderr, 'ERROR, if you want to add transcript info you must give both the --txalign AND --txsizes input'

# Create read objects
readDict = sizeDict(args.readsizes)

# Add adapter info to read objects 
with open(args.adaptalign, 'r') as f:
    # in the adapter blat file, the read is the target
    headFlag = False
    for line in f:
        # deal with header
        if line.startswith('psLayout'):
            headFlag = True
        elif line.startswith('------'):
            headFlag = False
            continue
        if headFlag:
            continue
        hit = BlatHit(line.strip())
        readObj = readDict[hit.tName] 
        # if match is in first half of the read, it's the left adapter
        if(hit.tSize/2 > hit.tEnd):
            readObj.addLA(hit.qName, hit.tEnd, hit.pslScore)
        # else it's the right adapter
        elif(hit.tSize/2 < hit.tStart):
            readObj.addRA(hit.qName, hit.tStart, hit.pslScore)
        # if the adapter crosses the middle of the read, there's something wrong
        else:
             print >>sys.stderr, "IGNORING: adapter {} for read {}, size {} maps to middle of read".format(hit.qName, hit.tName, hit.tSize)

if args.txalign:
    # Keep track of transcript sizes. If the input is in sam format, a fasta file with transcripts is required.
    txSizes = dict()
    # transcript hits can be sam or psl format
    if(args.txsizes):
        print >> sys.stderr, "assuming {} is in SAM format...".format(args.txalign)
        txSizes = txSizeDict(args.txsizes)
        with open(args.txalign, 'r') as f:
            for line in f:
                if line.startswith('@SQ') or line.startswith('@PG') or line.startswith('@HD'):
                    continue
                hit = SamHit(line.strip())
                readObj = readDict[hit.qName] 
                readObj.addTx(hit.tName, hit.qStart, hit.qEnd, hit.strand)
    
    # psl format
    else:
        # in the transcript file, one hit per read is expected and the read is the query
        print >> sys.stderr, "assuming {} is in PSL format...".format(args.txalign)
        with open(args.txalign, 'r') as f:
            headFlag = False
            for line in f:
                # deal with header
                if line.startswith('psLayout'):
                    headFlag = True
                elif line.startswith('------'):
                    headFlag = False
                    continue
                if headFlag:
                    continue
                hit = BlatHit(line.strip())
                txSizes[hit.tName] = hit.tSize 
                readObj = readDict[hit.qName] 
                readObj.addTx(hit.tName, hit.qStart, hit.qEnd, hit.strand)


# Print results
if args.txsizes:
    print "ReadID\tLeftAdaptEnd\tTxMatchStart\tTxMatchEnd\tRightAdaptStart\tReadSize\tleftAdapt\trightAdapt\tTranscript\tTxSize\tstrand"
else:
    print "ReadID\tLeftAdaptEnd\tRightAdaptStart\tReadSize\tleftAdapt\trightAdapt"
for r in readDict.values():
    r.printRead(args.txsizes)

