#
# A simple tool to generate reads with perfect coverage out a fasta stream.
#
# It creates paired end output from fragments of the same with half the reads
# coming from forward and half from reverse strands.
#
# It loads the whole genome into memory so might not work all that
# well for large genomes.
#
# Usage:
#
#    cat genome.fa | python perfect_coverage.py
#
# Will create the files called R1.fq and R2.fq
#
import sys, random

# How many reads per genomic coordinate.
COV = 1

# Read length
LEN = 50

# Fragment size
FRAG = 500

REVCOMP = dict(A='T', T='A', C='G', G='C')


def rc(x):
    return REVCOMP[x]


def revcomp(seq):
    out = ''.join(reversed(list(map(rc, seq))))
    return out


# Concatenate all lines of a FASTA file into
# a list of strings where each element represents
# a single full sequence.
genome, chrom = [], []
for line in sys.stdin:
    if line[0] == '>':
        if chrom:
            genome.append(''.join(chrom))
        chrom = []
    else:
        chrom.append(line.strip())

# Append last chromosome.
genome.append(''.join(chrom))

# Generate reads with the expected coverage.
read1, read2 = open('R1.fq', 'wt'), open('R2.fq', 'wt')

# Need unique names within a file.
pattern = "@P{}-{}\n{}\n+\n" + 'I' * LEN + '\n'

# Paired end reads
for chrom in genome:
    for index in range(len(chrom) - FRAG):

        # This is the first fragment.
        frag = chrom[index: index + FRAG]
        fwd_start, fwd_end = frag[:LEN], revcomp(frag[-LEN:])

        # Its reverse complement.
        rev_start, rev_end = revcomp(fwd_start), revcomp(fwd_end)

        # Generate the reads.
        # Turn indices into one based coordinates to read them more easily.
        for count in range(COV):
            read1.write(pattern.format(index + 1, count, fwd_start))
            read2.write(pattern.format(index + 1, count, fwd_end))
