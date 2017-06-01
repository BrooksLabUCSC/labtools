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

parser = argparse.ArgumentParser(description="From input sam file, get alignment stats such as fraction of query aligned and indels/mismatches from the CIGAR string and output histograms. Note: Reads must not be hard clipped.")

# Compulsory arguments start without a dash:
parser.add_argument('sam', type=str,  help="sam format file")
# optional boolean argument 
parser.add_argument('-q', '--quiet', help="Do not print warnings about skipped reads", action='store_true')
# optional regular argument
parser.add_argument('-b', '--outputbase', help="Output basename")

#TODO: Clarify % of alignment to % of aligned read
       
class SamHit(object):
    """
    Holds one sam format hit
    """
    def __init__(self, inline):
        fields = inline.split('\t')
        self.qname = fields[0]
        self.cigar = fields[5]
        # make sure all nucleotides are uppercase
        self.seq = fields[9].upper()
        self.charList = parseCigar(self.cigar)

    def stats(self):
        """Add information on the alignment"""

        # sanity check
        # at this point we should not have any hard clipped sequences
        self.size = len(self.seq)
        if self.size < 2:
            print >>sys.stderr, ("ERROR, read {} has no sequence, cannot calculate stats".format(self.qname))
            sys.exit(1)
        if self.cigar == '*':
            self.rawAlign = 0
            return

        if self.size != self.charList['M'] + self.charList['I'] + self.charList['S']:
            print >>sys.stderr,  self.qname, "ERROR: no match", self.size,  self.charList['M'], self.charList['I'], self.charList['S']
            sys.exit(1)

        self.nfract, self.gc =  ngcount(self.seq)

        # The raw alignment length is M+I
        self.rawAlign = self.charList['M'] + self.charList['I']

        # percent identity is matches divided by rawAlign
        self.pID = int(100*(float(self.charList['M']) / float(self.rawAlign)))

        # alignment percentage is rawAlign divided by self.size
        self.pAlign = int(100*(float(self.rawAlign)/float(self.size)))

        # counting indels as fraction of the length of the alignment or read is a bit treacherous, because we're counting
	# bases that were not part of the read (insertions) or not part of the genome (deletions)
        # Decision: Let's use the same rawAlign (matches plus insertions) as the divider. This makes it somewhat consistent with 
	# what we did before, even though the deletions are not part of this count
        self.ins = cigarFract(self.charList, 'I', self.rawAlign)
        self.delet = cigarFract(self.charList, 'D', self.rawAlign)
        self.indel = self.ins + self.delet 

def getBest(sams, args):
    """Return single hit with best match, replace sequence with full length if necessary"""
    if len(sams) == 1:
        return sams[0]
    # we're assuming no empty cigar strings if there are multiple alignments
    fullseq = False
    bestscore = 0
    for s in sams:
        if s.charList['H'] == 0:
            fullseq = s.seq
        else:  
            # put hard clip count into soft clip placeholder
            s.charList['S'] = s.charList['H']
        score = s.charList['M'] + s.charList['I']
        if score > bestscore:
            bestsam = s
            bestscore = score
    if not fullseq:
        if not args.quiet:
            print >>sys.stderr, ("WARNING, all sequences for {} appear to be hard clipped, skipping".format(s.qname))
        # modify the global variable
        global ignoredCount
        ignoredCount += 1
    bestsam.seq = fullseq
    return bestsam


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
    if cigar == '*':
        return charList
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


def yieldBestHit(f, sortsam, args):
    """Checks header for correct sort, then yields best hit for every read"""
    sortflag = False
    for line in f:
        if line.startswith('@HD'):
            if line.strip().split(':')[-1] == 'queryname':
                sortflag = True
            continue
        elif line.startswith('@'):
            continue
        break
    if not sortflag:
        # Instead of whining about it, just sort the file and start over
        # We're not using python's sort() function because the linux commandline sort is way more memory efficient
	# and our current file may be large.
        print >>sys.stderr, "No sort flag found in header, or not sorted by query. Sorting file to be sure." 
        with open(sortsam, 'w') as outf:
            for line in f:
                outf.write(line)
        f.close()
        # sort in place using the option -o
        sortcmd = 'sort {sfile} -o {sfile}'.format(sfile = sortsam)
        os.system(sortcmd)
        f = open(sortsam, 'r')
        # we have to read the first line in order to match the original file
        line = f.readline()
        phit = SamHit(line.strip())     
        readHits = [phit]
    else:
        # we're now in the body,  start yielding hits
        phit = SamHit(line.strip())     
        readHits = [phit]
        # Gather all alignments for a read
        # Note: bwa will not output the full sequence when another hit has it, so 
        # if there's a hard clip, we must retrieve the sequence, possibly from a later hit
    for line in f:
        hit = SamHit(line.strip())
        if hit.qname == phit.qname:
            readHits.append(hit)
        else:
            rhit = getBest(readHits, args)
            if rhit.seq:
              rhit.stats()
              yield rhit
            phit = hit
            readHits = [phit]
    rhit = getBest(readHits, args)
    if rhit.seq:
      rhit.stats()
      yield rhit


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
ignoredCount = 0

# We might need to sort the file: output it here
sortsam = ('.').join([outbase, str(os.getpid()), 'sorted.sam'])

# read the sam input file
with open(args.sam, 'r') as f:
    # skip header lines
    for hit in yieldBestHit(f, sortsam, args):
        hitCount += 1
        if hit.rawAlign == 0:
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

# Delete the sortfile if it's there
try:
    os.remove(sortsam)
except OSError:
    pass


# Now we have read the whole file so output some stats

uniqueReadCount = hitCount + ignoredCount
print '{0} of {1} reads are unaligned: {2:.1f}%'.format(unalignedCount, uniqueReadCount, 100* float(unalignedCount)/float(uniqueReadCount))
print '{0} of {1} reads were skipped: {2:.1f}%'.format(ignoredCount, uniqueReadCount, 100* float(ignoredCount)/float(uniqueReadCount))
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


