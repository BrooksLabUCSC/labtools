#!/usr/bin/python

import sys, os, re, argparse, textwrap, subprocess, shutil
sys.path.append('/pod/home/jeltje/lib')
from Bio import SeqIO
from collections import defaultdict
from intervaltree import Interval, IntervalTree

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\

This program corrects splices in (nanopore read) alignments (psl format), using a whole genome annotation and a junctions file (genepred format).
The junctions file should be derived from a short read alignment from the same sample, and can be created using
junctionsFromSam.py, which is based on 
https://raw.githubusercontent.com/anbrooks/juncBASE/master/preProcess_getASEventReadCounts.py

(NOTE TO SELF: junctions are reported but not used for any decisions yet. This could be optional.)
NOTE TO USER: If you don't have a junctions file, just use the genome annotation file as input for both the -j and -a parameters.

Before looking for splices, small gaps in the alignment are merged. The gap size can be changed; default is 30.
Splices are corrected if they are within a 'wiggle' distance of the intron start or intron end in the genome annotation.

For any novel splices, consensus sequences (GT/AG etc) are checked using the genome sequence. If no consensus is found, output junctions have '.' in the strand field.

The mRnaToGenes program must be in $PATH

OUTPUTS:
    novel.txt        contains a list of novel junctions with enough supporting reads
    junctions.bed    contains all junctions found in the query file. The score field contains the number of alignments with this junction.
    notfound.txt     is a list of query donor and acceptor sites that could not be found within the allowed wiggle distance
    multihit.txt     is a list of query donor and acceptor sites that had multiple hits within the allowed wiggle distance
    corrected.gp     contains all query annotations, splice corrected where possible

Notes: Novel junctions are only reported if they are identical in at least 3 annotations. This means that it is possible to miss junctions for which alignments are close but not identical. To see all novel junctions, set --novelthreshold to 1


        '''))
group = parser.add_argument_group('required arguments')
group.add_argument('-a','--annotations', type=str, required=True, help="genepred format genome annotation")
group.add_argument('-j','--junctions', type=str, required=True, help="genepred format junctions from SAM file")
group.add_argument('-g','--genofasta', type=str, required=True, help="genome fasta file")
group.add_argument('-q','--query', type=str, required=True, help="psl format alignment to be corrected")
group.add_argument('-w', '--wiggle', type=int, required = True, help="wiggle room for annotated splice junction")
group.add_argument('-o', '--outdir', type=str, default = '.', help="output directory")
parser.add_argument('-m', '--mergesize', type=int, default=30,  help="merge genome alignment gaps of this size (30)")
parser.add_argument('-n', '--novelthreshold', type=int, default=3,  help="report any novel junctions that are confirmed by at least this many reads")


class GpHit(object):
    """
    holds one genepred format alignment; extracts intron locations
    """
    def __init__(self, inline):
        fields = inline.split('\t')
        self.printBase = ('\t').join(fields[:8])
        self.qName, self.chrom, self.strand = map(str, fields[:3])
        self.genoStart, self.genoEnd, self.cdsStart, self.cdsEnd, self.blockCount = map(int, fields[3:8])
        self.chromStarts = fields[8].strip(',').split(',')
        self.chromEnds = fields[9].strip(',').split(',')
        self.lefts = map(int, self.chromEnds[:-1])
        self.rights = map(int, self.chromStarts[1:])
        self.makeIntrons() 

    def makeIntrons(self):
        self.introns = []
        for i in xrange(len(self.lefts)):
            intron = ('{}-{}').format(self.lefts[i], self.rights[i])
            self.introns.append(intron)

    def doPrint(self, ohandle):
        """
        Print current alignment in genePred format
        """
        if self.blockCount == 1:
            ohandle.write('{}\t{},\t{},\n'.format(self.printBase, self.chromStarts[0], self.chromEnds[0]))
            return
        chromStarts = str(self.genoStart)
        for i in self.rights:
            chromStarts += ',' + str(i)
        chromEnds = '' 
        for i in self.lefts:
            chromEnds += str(i) + ','
        chromEnds += str(self.genoEnd)
        ohandle.write('{}\t{},\t{},\n'.format(self.printBase, chromStarts, chromEnds))
        
def correctCoord(chrom, coord, jtree, wiggle, outnf, outmh):
    """
    Compare coord (integer) to intervals in jtree. If exact integer is found, return. If integer falls within range, 
    return middle of range. If integer falls within multiple ranges, return original integer and print warning.
    if integer is not found, return original and print warning

    """
    perfect, found = overlaps(jtree, coord, wiggle)
    if perfect:
        return coord
    if len(found) == 0:
       outnf.write("{}:{}\n".format(chrom, coord))
       return coord
    if len(found) >1:
       outmh.write("{}:{}\n".format(chrom, coord))
       return coord
    return found[0]


class Junctions():
    """
        Parse genepred format junctions 
    """
    def __init__(self):
        self.lefts = IntervalTree()
        self.rights = IntervalTree()
        self.introns = set()
        self.istrands = dict()
    def add(self, line, wiggle):
        """
           Add starts and ends as ranges in an intervaltree (faster for querying)
        """
        fields = line.strip().split('\t');

        chromStarts = fields[8].strip(',').split(',')
        rightset = map(int, chromStarts[1:])
        chromEnds = fields[9].strip(',').split(',')
        leftset = map(int, chromEnds[:-1])
        strand = fields[2]
        for i in xrange(len(leftset)):
            intron = ('{}-{}').format(leftset[i], rightset[i])
            self.iadd(strand, intron)
            self.lefts[leftset[i] - wiggle: leftset[i] + wiggle +1] = intron
            self.rights[rightset[i] - wiggle: rightset[i] + wiggle +1] = intron
            self.introns.add(intron)
    def iadd(self, strand, intron):
        """
        Keep track of intron associated strand
        """
        if intron in self.istrands:
            if self.istrands[intron] != strand:
                print >>sys.stderr, "WARNING, strand ambiguity for annotated junction", intron,  self.istrands[intron], strand
                self.istrands[intron] = '.'
        else:
            self.istrands[intron] = strand


def splitByChr(genepred, outdir):
    """
    Split input genepred file into one file per chromosome, named chrN.gp
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    chrList = []
    curChr = None
    foundchrs = set() 
    
    with open(genepred, 'r') as gp:
        for line in gp:
            chrom =line.split('\t')[1]
            if not curChr == chrom:
                if curChr:
                    with open(os.path.join(outdir, curChr + '.gp'), 'a') as o:
                        o.write(''.join(chrList))
                chrList = []
                curChr = chrom
                foundchrs.add(chrom)
            chrList.append(line)  
        with open(os.path.join(outdir, curChr + '.gp'), 'a') as o:
            o.write(''.join(chrList))
    return foundchrs

def overlaps(tree, coord, wiggle):
   """
    Returns the junction(s) that overlap
   """
   found = set()
   jcts = tree.search(coord)
   for i in jcts:
        if i.begin + wiggle == coord:
            return True, found
        found.add(i.begin + wiggle)
   return False, list(found)

def which(prog):
    """See if a program exists on the system and return the path"""
    cmd = ["which",prog]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    res = p.stdout.readline().rstrip()
    if len(res) == 0: 
        print >>sys.stderr, 'ERROR, cannot find program {}, please install and/or add to $PATH'.format(prog)
        sys.exit(1)
    return res

def listCount(mylist, wiggle, cutoff):
    """
    Count occurrences of coordinate in list within wiggle distance, print if more than cutoff occurrences
    """
    mylist = sorted(mylist)
    cur = 0
    hits = []
    for i in mylist:
        if cur + wiggle >= i:
            hits.append(i)
        else:
            if len(hits) >= cutoff:
                print >>sys.stderr, hits
            hits = []
        cur = i

def makeBed(chr, chr_seq, intron, score, outfile, strandlist):
    """
    Create bed format intron output with a blocksize of 20
    If the intron is not found in existing annotation, or if strand info conflicts, 
    infers strand from intron start and end sequence 
    """
    start, end = map(int, intron.split('-'))
    strand = '.'
    if intron in strandlist:
        strand = strandlist[intron]
    if strand == '.':
        strand = disambiguateJcnStr(chr_seq, start, end)

    cstart = start -20
    cend = end +20
    bedline = '{chr}\t{cstart}\t{cend}\t{chr}:{intron}\t{score}\t{strand}\t{cstart}\t{cend}\t0,0,255\t2\t20,20\t0,{bstart},\n'.format(chr=chr, 
         strand=strand, cstart=cstart, cend=cend, intron=intron, score=score, bstart=end-cstart)
    outfile.write(bedline)


def disambiguateJcnStr(chr_seq, start, end):
    """
    Will use splice site sequence to infer strand
    If no strand can be determined, returns '.'
    This function is from  disambiguate_junctions.py by Angela Brooks
    """

    strand = '.'

    # extract the sequence from the chromosome
    intron_seq = chr_seq[start-1:end]

    if intron_seq.startswith("GT") and intron_seq.endswith("AG"):
        strand = "+"
    elif intron_seq.startswith("CT") and intron_seq.endswith("AC"):
        strand = "-"
    # Other common splice site sequence
    elif intron_seq.startswith("GC") and intron_seq.endswith("AG"):
        strand = "+"
    elif intron_seq.startswith("CT") and intron_seq.endswith("GC"):
        strand = "-"
    # minor spliceosome
    elif intron_seq.startswith("AT") and intron_seq.endswith("AC"):
        strand = "+"
    elif intron_seq.startswith("GT") and intron_seq.endswith("AT"):
        strand = "-"
    # Priority to 5' splice site since there is more information
    # there
    elif intron_seq.startswith("GT"):
        strand = "+"
    elif intron_seq.endswith("AC"):
        strand = "-"
    elif intron_seq.endswith("AG"):
        strand = "+"
    elif intron_seq.startswith("CT"):
        strand = "-"
    else:
        print >>sys.stderr, "Cannot find strand for {}-{}".format(start, end)
    return strand


def findChr(records, this_chr):
    """
    Will identify chr sequence, then call disambiguateJcnStr
    This function is from  disambiguate_junctions.py by Angela Brooks
    """
    chr_seq = None
    if this_chr in records:
        chr_seq = records[this_chr].seq
    elif this_chr.lstrip("chr") in records:
        chr_seq = records[unFormatChr(this_chr)].seq
    else:
        print "Cannot find %s in genome sequence." % this_chr
        sys.exit(1)
    return chr_seq


# Main
args = parser.parse_args()

# Sanity check
if args.mergesize <= args.wiggle:
    print >>sys.stderr, "ERROR: mergesize (currently {}) must be larger than wiggle (currently: {}), please correct and rerun".format(args.mergesize, args.wiggle)
    sys.exit(1)

# Setup temporary output location
tmpdir = '/tmp/spliceCorr_' + str(os.getpid())
if os.path.exists(tmpdir):
    shutil.rmtree(tmpdir)
os.makedirs(tmpdir)

# Convert the genome alignment psl to genepred - doing this first because it is most likely to fail
mrnaToGene = which('mrnaToGene')
msize = '-insertMergeSize={}'.format(args.mergesize)
gpgeno = '{}/{}.gp'.format(tmpdir, args.query)
cmd = [mrnaToGene, '-noCds', msize, args.query, gpgeno]
try:
    subprocess.check_call(cmd)
except subprocess.CalledProcessError:
    print >>sys.stderr, "Cannot run mrnaToGene correctly, exiting."
    sys.exit(1)

# split all input files into chromosomes
workdir = args.outdir
if not os.path.exists(workdir):
    os.makedirs(workdir)
jdir = os.path.join(tmpdir, 'junctions')
adir = os.path.join(tmpdir, 'annots')
qdir = os.path.join(tmpdir, 'aligns')
splitByChr(args.junctions, jdir)
splitByChr(args.annotations, adir)
chroms = splitByChr(gpgeno, qdir)

# read genome fasta
try:
    records = SeqIO.index(args.genofasta, "fasta")
except:
    print >>sys.stderr, "ERROR: Could not open genome sequence."
    sys.exit(1)

# Per-chromosome analysis
with open(os.path.join(workdir, 'corrected.gp'), 'w') as outgp, \
  open(os.path.join(workdir, 'junctions.bed'), 'w') as outbed, \
  open(os.path.join(workdir, 'novel.txt'), 'w') as outnv, \
  open(os.path.join(workdir, 'notfound.txt'), 'w') as outnf, \
  open(os.path.join(workdir, 'multihit.txt'), 'w') as outmh:
    outnv.write("junction\talignment\tRNASeq\tannotation\n")
    for c in chroms:
        print >>sys.stderr, c
        cors = Junctions()
        junctionIntrons = []
        annotIntrons = []
        # read in junctions file
        with open(os.path.join(jdir, c + '.gp'), 'r') as f:
            for line in f:
                cors.add(line, args.wiggle)
                hit = GpHit(line.strip())
                junctionIntrons.extend(hit.introns)
        # add the annotation file
        with open(os.path.join(adir, c + '.gp'), 'r') as f:
            for line in f:
                cors.add(line, args.wiggle)
                hit = GpHit(line.strip())
                annotIntrons.extend(hit.introns)
        # read alignment file and correct intronstarts and intron ends
        alignIntrons = []
        with open(os.path.join(qdir, c + '.gp'), 'r') as f:
            for line in f:
                hit = GpHit(line.strip())
                newlefts = []
                for i in hit.lefts:
                    newlefts.append(correctCoord(c, i, cors.lefts, args.wiggle, outnf, outmh))
                newrights = []
                for i in hit.rights:
                    newrights.append(correctCoord(c, i, cors.rights, args.wiggle, outnf, outmh))
                hit.lefts = newlefts
                hit.rights = newrights
                hit.makeIntrons()
                alignIntrons.extend(hit.introns)
                hit.doPrint(outgp)
        
        chr_seq = findChr(records, c) 
        # now for every unique junction in the alignment, count the number of supporting reads and check if the junction is present in the annotation and in the junctions file
        # Also print out every junction in bed format and check for consensus splice sites.
        for i in (set(alignIntrons)):
            alcount = alignIntrons.count(i) 
            jucount = junctionIntrons.count(i)
            ancount = annotIntrons.count(i)
            makeBed(c, chr_seq, i, alcount, outbed, cors.istrands)
            if not ancount and alcount >= args.novelthreshold:
               outnv.write("{}:{}\t{}\t{}\t{}\n".format(c, i, alignIntrons.count(i), junctionIntrons.count(i), annotIntrons.count(i)))

shutil.rmtree(tmpdir)

print textwrap.dedent('''\

SUCCESS

These output files were created:
    novel.txt        contains a list of novel junctions with enough supporting reads
    junctions.bed    contains all junctions found in the query file. The score field contains the number of alignments with this junction.
    notfound.txt     is a list of query donor and acceptor sites that could not be found within the allowed wiggle distance
    multihit.txt     is a list of query donor and acceptor sites that had multiple hits within the allowed wiggle distance
    corrected.gp     contains all query annotations, splice corrected where possible

Notes: 
   Novel junctions are only reported if they are identical in at least 3 annotations. This means that it is possible to miss junctions for which alignments are close but not identical. To see all novel junctions, set --novelthreshold to 1
   If the bed file contains a period in the strand field, the junction is likely an artifact (for example a misalignment of exon ends or an unfilled gap in an exon)

PROGRAM FINISHED

        ''')

