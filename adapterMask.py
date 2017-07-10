#!/usr/bin/python

import sys, os, re, getopt, argparse, textwrap

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\

This program is based on nanoporeMatchTable py.
It takes in a psl file of adapter queries to a nanopore read database and outputs a bed file with coordinates
that can be used to mask the sequences.
Suggested blat settings:
	blat -minIdentity=0 -minMatch=0 -minScore=0

Matches within a settable limit from each end of the read (default 100) result in masking out the full end.
Matches outside this region with a set coverage (default 50nt of the read) are also masked out

        '''))
group = parser.add_argument_group('required arguments')
group.add_argument('--adaptalign', type=str, required=True, help="psl format alignment of ADAPTERS vs READS")
#group.add_argument('--readsizes', type=str, required=True, help="faCount output for nanopore reads")
# optional flags
#parser.add_argument('--minEndMatch', default=100, help="", action='store_true')
parser.add_argument('--minEndMatch', type=int, default=100, help="If the adapter matches within this range of the 5' or 3' end, the whole end is masked (default 100)")
parser.add_argument('--minCover', type=int, default=50, help="If the adapter matches away from the end, it must cover at least this much of the read to be output (default 50)")


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
        self.tSize, self.tStart, self.tEnd = map(int, fields[14:17])
        self.blockCount = int(fields[17])
        self.blockSizes, self.qStarts, self.tStarts = fields[18:]
        self.qCover = self.qEnd - self.qStart	# length of blat match
        self.pslScore = self.match + self.rep - self.mismatch - self.qGapCount - self.tGapCount #- self.tGapBases
        self.inline = inline
    def doPrint(self, endsize, cover):
        if hit.tEnd < endsize:
            print '{}\t{}\t{}'.format(hit.tName, 0, hit.tEnd)
        elif hit.tSize - hit.tEnd < endsize:
            print '{}\t{}\t{}'.format(hit.tName, hit.tEnd, hit.tSize)
        elif hit.tEnd - hit.tStart > cover:
            print '{}\t{}\t{}'.format(hit.tName, hit.tStart, hit.tEnd)


args = parser.parse_args()
# Run program

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
        hit.doPrint(args.minEndMatch, args.minCover)


