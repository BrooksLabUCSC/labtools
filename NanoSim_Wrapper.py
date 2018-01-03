#!/usr/bin/env python

import sys, os, glob, re, argparse, textwrap, subprocess, random
from collections import defaultdict
from natsort import natsorted

parser = argparse.ArgumentParser(formatter_class= argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''
    This program creates individual fasta files that will be used as input for NanoSim as well as a text file 
    that designates how many reads are wanted for each transcript along with the minimum and maximum values 
    for NanoSim commands '''))

group = parser.add_argument_group('required arguments')
group.add_argument('-c', '--count', type=str, required=True, help="Tab delimited file that contains desired transcript headers along with the number of desired reads")
group.add_argument('-f', '--fasta_file', type=str, required=True, help="fasta file that contains all the potential transcripts that could be simulated (gencode)")
group.add_argument('-p', '--mincov', type=float, default=0.9, help="Minimum transcript coverage to be considered full length (default 0.9)")
group.add_argument('-n', '--not_full_length', type=float, default=0.2, help="What fraction of the output reads should be not full length (default 0.2)")
group.add_argument('-M', '--model_prefix', type=str, required=True, help="Model prefix used to simulate reads e.g. Hg38_Model2/Hg38_training")
group.add_argument('-d', '--outdir', type=str, default="ns_output", help="Output directory")
group.add_argument('-O', '--combined', type=str, default="combined.fasta", help="final combined simulated reads fasta file")
#optional flag
group.add_argument('--min_len', type=int, default=0, help="Transcripts less than this length will not be simulated")
group.add_argument('--max_len', type=int, required=False, help="Transcripts greater than this length will not be simulated")
group.add_argument('-m', '--multiply', type=int, default=1, help="Use to increase the number of reads you want per transcript, will use the count number from count file and will multiply by number given. Use 1 to use original count number")

def getCounts(infile):

    counts = dict()
    with open(infile, 'r') as f:
        for line in f:
            line = line.strip()
            [id, score] = line.split("\t")
            counts[id] = int(score)
    return counts

def yieldFasta(fasta):
    """Generator function to yield fasta id and sequence"""
    ok = True
    try:
        with open(fasta, 'r') as fh:
            while ok:
                try:
                    line = next(fh)
                except:
                    ok = False
                    yield id, sequence
                if line[0] == '>':
                    try:
                        yield id, sequence
                    except NameError:
                        pass
                    # remove newline and '>'
                    id = line.strip('\n')
                    id = id[1:]
                    sequence = ''
                elif line[0] != '>':
                    sequence += line.strip('\n')

    except (IOError, OSError):
        print("Error opening / processing file")
    except StopIteration:
        pass

def runNanosimPerTx(args, countDict):
    """For every transcript in the input fasta, see if it's in the counts file and if it passes other requirements. If so, run NanoSim"""
    prog = which('simulator.py')
    for id, sequence in yieldFasta(args.fasta_file):
	if id in countDict:
            print 'simulating for {}'.format(id)
            seq_len = len(sequence)
            if seq_len < args.min_len:
                continue
            if args.max_len:
                if seq_len > args.max_len:
                    continue
            count = countDict[id] * args.multiply
            minimum = int(float(seq_len) * args.mincov)
            num_full_length = int(float(count) * (1- args.not_full_length))
            num_non_full_length = int(float(count) * args.not_full_length)
            # create inputfile for Nanosim
            nsinput = os.path.join(args.outdir, "tmp.fa")
            with open(nsinput, 'w') as output:
                output.write(">{}\n{}\n".format(id, sequence))
            # nanoSim creates two files, capture their names
            runNanoSim(prog, args, nsinput, num_full_length, num_non_full_length, seq_len, minimum)
        else:
            print 'skipping {}'.format(id)

def yieldCount():
    """Simple generator function to yield an incremental number"""
    i=0
    while True:
        i += 1
	yield i
        
def runNanoSim(NanoSim, args, inputfa, num_reads_full, num_reads_part, seq_size, min_len):

    # full length reads have a size range from min_len to seq_size
    outfileFull = "{}/full.{}".format(args.outdir, counter.next())
    cmd=[
        NanoSim, 'linear',
        '-r', inputfa,
        '-c', args.model_prefix,
        '-n', str(num_reads_full),
        '--max_len', str(seq_size),
        '--min_len', str(min_len),
        '-o', outfileFull
        ]        
    subprocess.call(cmd)

    # partial reads have a size range from 100 to min_len
    outfilePart = "{}/part.{}".format(args.outdir, counter.next())
    cmd=[
        NanoSim, 'linear',
        '-r', inputfa,
        '-c', args.model_prefix,
        '-n', str(num_reads_part),
        '--max_len', str(min_len -1),
        '--min_len', "100",
        '-o', outfilePart
        ]        
    subprocess.call(cmd)
    return [outfilePart, outfileFull]

def combineFa(args, read_files):

    with open(args.combined, "w") as outfile:
        for f in read_files:
            with open(f, "r") as infile:
                outfile.writelines(infile.read())



def which(prog):

    cmd = ["which", prog]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    res = p.stdout.readline().rstrip()
    if len(res) == 0:
        print >>sys.stderr, 'ERROR, cannot find program {}, please add to $PATH'.format(prog)
        sys.exit(1)
    return(res)


def main():
    args = parser.parse_args()
    if os.path.exists(args.outdir):
        print >>sys.stderr, 'ERROR, output directory {} already exists, please remove or rename'.format(args.outdir) 
        sys.exit(1)
    os.makedirs(args.outdir)

    # get counts per transcript
    countsDict = getCounts(args.count)
    # set a counter
    global counter 
    counter = yieldCount()
    # run NanoSim
    runNanosimPerTx(args, countsDict)
    outfiles = glob.glob("{}/*.fasta".format(args.outdir))
    # combine output
    combineFa(args, outfiles)

if __name__ == '__main__':
        main()

