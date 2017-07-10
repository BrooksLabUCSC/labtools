#!/usr/bin/python

import sys, os, re, argparse, textwrap, subprocess, shutil
#sys.path.append('/pod/home/jeltje/lib')
from collections import defaultdict
from intervaltree import Interval, IntervalTree

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\

Create bed output format junctions file

        '''))
group = parser.add_argument_group('required arguments')
group.add_argument('annotations', type=str, help="genepred format genome annotation")
#group.add_argument('-g','--annotations', type=str, required=True, help="genepred format genome annotation")
group.add_argument('-o','--outdir', type=str, default='.', help="temporary output dir")


class Junctions():
    """
        Parse genepred format junctions 
    """
    def __init__(self, line):
        self.introns = set()
#    def add(self, line):
        fields = line.strip().split('\t');

        chromStarts = fields[8].strip(',').split(',')
        rightset = map(int, chromStarts[1:])
        chromEnds = fields[9].strip(',').split(',')
        leftset = map(int, chromEnds[:-1])
        strand = fields[2]
        for i in xrange(len(leftset)):
            intron = ('{}_{}_{}').format(leftset[i], rightset[i], strand)
            self.introns.add(intron)

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


def makeBed(chr, intron, score):
    """
    Create bed format output with a blocksize of 20
    """
    start, end, strand = intron.split('_')
    start = int(start)
    end = int(end)
    cstart = start -20
    cend = end +20
    print '{chr}\t{cstart}\t{cend}\t{chr}:{start}-{end}\t{score}\t{strand}\t{cstart}\t{cend}\t0,0,255\t2\t20,20\t0,{bstart},'.format(chr=chr, 
         start=start+1, end=end, cstart=cstart, cend=cend, strand=strand, intron=intron, score=score, bstart=end-cstart)


# Main
args = parser.parse_args()

tmpdir = '/tmp/makesplice_' + str(os.getpid())
# split input file into chromosomes
workdir = args.outdir
if not os.path.exists(workdir):
    os.makedirs(workdir)
qdir = os.path.join(tmpdir, 'aligns')
chroms = splitByChr(args.annotations, qdir)

for c in chroms:
        #cors = Junctions()
        junctionIntrons = []
        # read in junctions file
        with open(os.path.join(qdir, c + '.gp'), 'r') as f:
            for line in f:
                #cors.add(line)
                hit = Junctions(line)
                junctionIntrons.extend(hit.introns)
        for i in (set(junctionIntrons)):
            makeBed(c, i, junctionIntrons.count(i))

shutil.rmtree(tmpdir)

