#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys, os, re, getopt, argparse
import subprocess
import numpy

# Argparse deals with getting input arguments from the user. 
# It also results in a short usage statement when the program is called with '-h'
# and a longer statement when called with '--help'

parser = argparse.ArgumentParser(description="From input sam file, get alignment stats such as fraction of query aligned and indels/mismatches from the CIGAR string and output histograms.")

# Compulsory arguments start without a dash:
parser.add_argument('sam', type=str,  help="sam format file")
# optional boolean argument (not used here)
# parser.add_argument('-d', '--debug', help="Optional debugging output", action='store_true')
# optional regular argument
parser.add_argument('-b', '--outputbase', help="Output basename")

#TODO: Allow bam input; Sort input on read ID (samtools does this); take one alignment only
#TODO: Clarify % of alignment to % of aligned read
#TODO: Allow hard clip input (cannot calculate GC% if used) 
       
class SamHit(object):
    """
    Holds one sam format hit
    """
    def __init__(self, inline):
        fields = inline.split('\t')
        cigar = fields[5]
### RETHINK THIS - doesn't work with all sam inputs #####
        # make sure all nucleotides are uppercase
        seq = fields[9].upper()
        self.size = len(seq)
        if self.size == 0:
            self.unaligned = True
            return
        self.nfract, self.gc =  ngcount(seq)
    
###########################
        if cigar == '*':
            self.unaligned = True
            return
        self.unaligned = False
        charList = parseCigar(cigar)
        # sanity check - this fails on the hard clipped reads. We should not have hard clips.
        if self.size != charList['M'] + charList['I'] + charList['S'] + charList['H']:
            print >>sys.stderr,  fields[0], "no match", self.size,  charList['M'], charList['I'], charList['S'], charList['H']

        # The raw alignment length is M+I
        self.rawAlign = charList['M'] + charList['I']

        # percent identity is matches divided by rawAlign
        self.pID = int(100*(float(charList['M']) / float(self.rawAlign)))

        # alignment percentage is rawAlign divided by self.size
        self.pAlign = int(100*(float(self.rawAlign)/float(self.size)))

        # indels are counted a bit differently, we use the length of the read minus introns and overhang
        alignSize = self.size - (charList['S'] + charList['N'] + charList['H'])
        
        self.ins = cigarFract(charList, 'I', alignSize)
        self.delet = cigarFract(charList, 'D', alignSize)
        self.indel = self.ins + self.delet 

def ngcount(seq):
    """Returns N fraction and GC content of input sequence"""
    ncount = seq.count('N')
    if ncount == len(seq):
        nfract = 100
        gcfract = 0
    else:
        nfract = float(ncount)/float(len(seq))
        eff = float(len(seq) - ncount)
        gc = float(seq.count('G') + seq.count('C'))
        gcfract = int(100 * gc/eff)
    return nfract, gcfract

def parseCigar(cigar):
    """Return a dictionary with cigar characters and their total length"""
    # M is match, D is deletion of reference base in read, I is insertion, N is intron (equivalent to D), S is overhang
    # H is hard clip

    # First identify all CIGAR letters. We know these exist:
    charList = { 'M':0, 'I':0, 'D':0, 'S':0, 'N':0, 'H':0 } 
    # This gets all letters in the current cigar
    chars = list(int_filter(set(list(cigar))))
    for c in chars:
        charList[c] = 0

    # Now count how many total matches, mismatches, etc
    searchString = "({})([{}])".format('\d+', ('').join(c for c in chars))
    matches = re.findall(searchString, cigar)
    # and add total counts for each
    for m in matches:
        charList[m[1]] += int(m[0])
    return charList


def cigarFract(charList, char, size):
    """Return fraction of given symbol in total cigar string"""
    if charList[char] > 0:
        return int(100 * (float(charList[char])/float(size)))
    else:
        return 0
            


def int_filter(myList):
    """Return value if it's not an integer, can be used to remove numbers from CIGAR"""
    for v in myList:
        try:
            int(v)
            continue # Skip these
        except ValueError:
            yield v # Keep these

def sizePlot(mylist, color, label, outfile):
    """Create histogram of read sizes"""
    plt.hist(mylist, bins=range(0, max(mylist), 100), facecolor=color, label=label, alpha=0.75)

    if outfile:
        plt.xlabel('Size(nt)')
        plt.ylabel('Count') 
        title = 'Read and alignment lengths' #'{} reads'.format(label, len(mylist))
        plt.title(title)
        plt.legend()
        plt.grid(True)
        plt.savefig(outfile)     

def indelPlot(mylist, color, label, outfile):
    """Create histogram with insertions and deletions"""
    val, weight = zip(*[(k, v) for k,v in mylist.items()])
    plt.hist(val, weights=weight, bins=range(0, max(mylist), 1), facecolor=color, alpha=0.75, label=label)
    if outfile:
        plt.xlabel('% of read covered')
        plt.ylabel('Count') 
        title = 'Read indels' 
        plt.title(title)
        plt.grid(True)
        plt.legend()
        plt.savefig(outfile)     

def fractPlot(mylist, color, label, outfile):
    """Create histogram (appends to same plot when called multiple times)"""
    val, weight = zip(*[(k, v) for k,v in mylist.items()])
    plt.hist(val, weights=weight, bins=range(0, 101, 1), facecolor=color, alpha=0.75, label=label)
    if outfile:
        plt.xlabel('% of read covered')
        plt.ylabel('Count') 
        title = 'Read QC' #'{} reads'.format(label, len(mylist))
        plt.title(title)
        plt.grid(True)
        plt.legend(loc='upper center')
        plt.savefig(outfile)     

# Main
# read in command line and options

args = parser.parse_args()
outbase = args.outputbase or './samstats'
# Run program

# counters for each value (0-100) of each feature - takes less memory than storing in lists.
gcList = dict()
insList = dict()
deletList = dict()
indelList = dict()
palignList = dict()
pidList = dict()
sizeList = []
rawalignList = []
hitCount = 0
unalignedCount = 0

# read the sam input file
with open(args.sam, 'r') as f:
    # skip header lines
    for line in f:
        if line.startswith('@'):
            continue
        # Object call to samHit gets all the important info for each read
        hit = SamHit(line.strip())
        hitCount += 1
        if hit.unaligned == True:
            unalignedCount += 1
            continue
        # add stats for the current read to our counters. All fractions are in whole percentage points.
        gcList[hit.gc] = gcList.get(hit.gc,0) +1
        insList[hit.ins] = insList.get(hit.ins,0) +1
        deletList[hit.delet] = deletList.get(hit.delet,0) +1
        indelList[hit.indel] = indelList.get(hit.indel,0) +1
        palignList[hit.pAlign] = palignList.get(hit.pAlign,0) +1
        pidList[hit.pID] = pidList.get(hit.pID,0) +1
        sizeList.append(hit.size)
        rawalignList.append(hit.rawAlign)

# Now we have read the whole file so output some stats
print '{} of {} reads are unaligned: {}%'.format(unalignedCount, hitCount, unalignedCount/hitCount)
print >>sys.stderr, "Done counting, creating figures..."

# And create the histograms
outfile = ('.').join([outbase, 'sizes.png'])
print >>sys.stderr, outfile
sizePlot(rawalignList, 'green', 'alignment lengths', False)
sizePlot(sizeList, 'red', 'read sizes', outfile)

# After we print the first plot, clean the plt object
plt.clf()
plt.cla()

# percentages
outfile = ('.').join([outbase, 'QC.png'])
print >>sys.stderr, outfile
fractPlot(gcList, 'green', 'GC', False)
fractPlot(indelList, 'blue', 'indels', False)
fractPlot(palignList, 'purple', 'percent aligned', False)
fractPlot(pidList, 'orange', 'percent identity', outfile)

plt.clf()
plt.cla()

# indels
outfile = ('.').join([outbase, 'indel.png'])
print >>sys.stderr, outfile
indelPlot(deletList, 'blue', 'deletions', False)
indelPlot(insList, 'yellow', 'insertions', outfile)


