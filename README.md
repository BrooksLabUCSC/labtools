# labtools
Small standalone scripts


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
