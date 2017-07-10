# labtools
Small standalone scripts

**correctSplice.py**
Corrects splices in (nanopore read) alignments (psl format), using a whole genome annotation and a junctions file (genepred format).

**nanoporeQC**
Bash program that runs blat to align adapters to nanopore reads and formats input for `nanoporeMatchTable.py`, then calls the program

**nanoporeMatchTable.py**
Creates a table with adapter positions in nanopore reads, based on a `psl` format input file.

**samAlignStats.py**
Creates a set of pretty histograms to show information on alignments, such as percentage indels and mismatches.

**adapterMask.py**
Creates a bed file for masking adapters, which can be used as input to [maskOutFa](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/maskOutFa)

**genepredToJunctions.py**
Read a genepred input file, output all junctions in bed format

**junctionsFromSam.py**
Takes a SAM format file as input and creates an output directory with a junction file in bed format.

# Details on labtools

## correctSplice.py
This program corrects splices in (nanopore read) alignments (psl format), using a whole genome annotation and a junctions file (genepred format).
```
usage: correctSplice.py [-h] -a ANNOTATIONS -j JUNCTIONS -g GENOFASTA -q QUERY
                        -w WIGGLE [-o OUTDIR] [-m MERGESIZE]
                        [-n NOVELTHRESHOLD]
optional arguments:
  -h, --help            show this help message and exit
  -m MERGESIZE, --mergesize MERGESIZE
                        merge genome alignment gaps of this size (30)
  -n NOVELTHRESHOLD, --novelthreshold NOVELTHRESHOLD
                        report any novel junctions that are confirmed by at
                        least this many reads

required arguments:
  -a ANNOTATIONS, --annotations ANNOTATIONS
                        genepred format genome annotation
  -j JUNCTIONS, --junctions JUNCTIONS
                        genepred format junctions from SAM file
  -g GENOFASTA, --genofasta GENOFASTA
                        genome fasta file
  -q QUERY, --query QUERY
                        psl format alignment to be corrected
  -w WIGGLE, --wiggle WIGGLE
                        wiggle room for annotated splice junction
  -o OUTDIR, --outdir OUTDIR
                        output directory
```
The junctions file should be derived from a short read alignment from the same sample, and can be created using
junctionsFromSam.py, which is based on 
https://raw.githubusercontent.com/anbrooks/juncBASE/master/preProcess_getASEventReadCounts.py

NOTE TO USER: If you don't have a junctions file, just use the genome annotation file as input for both the -j and -a parameters.

Before looking for splices, small gaps in the alignment are merged. The gap size can be changed; default is 30.
Splices are corrected if they are within a 'wiggle' distance of the intron start or intron end in the genome annotation.

For any novel splices, consensus sequences (GT/AG etc) are checked using the genome sequence. If no consensus is found, output junctions have '.' in the strand field.

**The mRnaToGenes program must be in $PATH**

```
OUTPUTS:
    novel.txt        contains a list of novel junctions with enough supporting reads
    junctions.bed    contains all junctions found in the query file. The score field contains the number of alignments with this junction.
    notfound.txt     is a list of query donor and acceptor sites that could not be found within the allowed wiggle distance
    multihit.txt     is a list of query donor and acceptor sites that had multiple hits within the allowed wiggle distance
    corrected.gp     contains all query annotations, splice corrected where possible
```
Notes: Novel junctions are only reported if they are identical in at least 3 annotations. This means that it is possible to miss junctions for which alignments are close but not identical. To see all novel junctions, set --novelthreshold to 1

## nanoporeQC
This bash program runs blat to align adapters to nanopore reads and formats input for `nanoporeMatchTable.py`, then calls the program

Prerequisites are [blat](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat) and [faCount](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faCount)

Inputs are a reads file (`fasta` or `fastq`), and an adapter file (`fasta` with two sequences)

When ALSO given a reads-vs-transcripts alignment file (in `sam` or `psl` format), 
it adds the transcript matches to the output table

IF you add a sam format input, you MUST ALSO add a fasta file with the transcripts
that were used in the alignments

The inputs must be in order, like so:

```
nanoporeQC reads.fa adapters.fa > output.tsv
nanoporeQC reads.fq adapters.fa alignments.psl > output.tsv
nanoporeQC reads.fa adapters.fa alignments.sam transcripts.fa > output.tsv
```

## nanoporeMatchTable.py
This program creates a table with adapter positions in nanopore reads, based on a `psl` format input file.
Optionally a second alignment file of reads vs transcripts can be added in `psl` or `sam` format, the program
then adds start and end positions of the transcript alignment in the read.
If the `sam` format is used, you MUST give in a file with transcript sizes in the format

`TranscriptID	size`

The program [faCount](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faCount) outputs this format (you can leave the header and footer).

```
optional arguments:
  -h, --help            show this help message and exit
  --txalign TXALIGN     psl or sam format file of READS vs TRANSCRIPTS
  --txsizes TXSIZES     faCount output for transcripts, required if using sam
                        format

required arguments:
  --adaptalign ADAPTALIGN
                        psl format alignment of ADAPTERS vs READS
  --readsizes READSIZES
                        faCount output for nanopore reads
```
## samAlignStats.py
This program creates a set of pretty histograms to show information on alignments, such as percentage indels and mismatches.

```
usage: samAlignStats.py [-h] [-q] [-b OUTPUTBASE] sam

From input sam file, get alignment stats such as fraction of query aligned and
indels/mismatches from the CIGAR string and output histograms. Note: Reads
must not be hard clipped.

positional arguments:
  sam                   sam format file

optional arguments:
  -h, --help            show this help message and exit
  -q, --quiet           Do not print warnings about skipped reads
  -b OUTPUTBASE, --outputbase OUTPUTBASE
                        Output basename
```

## adapterMask.py
Creates a bed file for masking adapters, which can be used as input to [maskOutFa](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/maskOutFa)

```
usage: adapterMask.py [-h] --adaptalign ADAPTALIGN [--minEndMatch MINENDMATCH]
                      [--minCover MINCOVER]

This program is based on nanoporeMatchTable py.
It takes in a psl file of adapter queries to a nanopore read database and outputs a bed file with coordinates
that can be used to mask the sequences.
Suggested blat settings:
	blat -minIdentity=0 -minMatch=0 -minScore=0

Matches within a settable limit from each end of the read (default 100) result in masking out the full end.
Matches outside this region with a set coverage (default 50nt of the read) are also masked out

optional arguments:
  -h, --help            show this help message and exit
  --minEndMatch MINENDMATCH
                        If the adapter matches within this range of the 5' or
                        3' end, the whole end is masked (default 100)
  --minCover MINCOVER   If the adapter matches away from the end, it must
                        cover at least this much of the read to be output
                        (default 50)

required arguments:
  --adaptalign ADAPTALIGN
                        psl format alignment of ADAPTERS vs READS
```

## genepredToJunctions.py
Read a genepred input file, output all junctions in bed format

```
usage: genePredToJunctions.py [-h] [-o OUTDIR] annotations

Create bed output format junctions file

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  annotations           genepred format genome annotation
  -o OUTDIR, --outdir OUTDIR
                        temporary output dir
```

## junctionsFromSam.py
Takes a SAM format file as input and creates an output directory with a junction file in bed format.
This program is an adaptation of preProcess_getASEventReadCounts.py in [juncBase](https://github.com/anbrooks/juncBASE.git)
```
Usage: junctionsFromSam.py [options]

Options:
  -h, --help           show this help message and exit
  -s SAM_FILE          SAM/BAM file of read alignments to junctions and
                       the genome. More than one file can be listed,
                       but comma-delimited, e.g file_1.bam,file_2.bam
  --unique             Only keeps uniquely aligned reads. Looks at NH
                       tag to be 1 for this information.
  -n NAME              Name prefixed to output files and used for
                       output BED file. Default=getASEventReadCounts_input
  -l READ_LENGTH       Expected read length if all reads should be of
                       the same length
  -c CONFIDENCE_SCORE  The mininmum entropy score a junction
                       has to have in order to be considered
                       confident. The entropy score =
                       -Shannon Entropy. Default=1.0
  -j FORCED_JUNCTIONS  File containing intron coordinates
                       that correspond to junctions that will be
                       kept regardless of the confidence score.
  -o OUTPUT_DIR        Directory to place all output files.

```
