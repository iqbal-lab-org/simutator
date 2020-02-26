[![Build Status](https://travis-ci.com/iqbal-lab-org/simutator.svg?branch=master)](https://travis-ci.com/iqbal-lab-org/simutator)

# simutator

Simulate mutations in genomes.

## Dependencies

Python3 (known to work with version 3.6.9), and optionally the read simulator
[ART](http://bioinformatics.oxfordjournals.org/content/28/4/593.abstract)
(only needed if you simulated reads, not needed to mutate genomes).
To simulate reads, the program `art_illumina` must be in your `$PATH`.
If you use ART, please cite "ART: a next-generation sequencing read simulator",
Huang et al 2011.

## Install

Git clone this repository and run:

```
    pip3 install .
```

It will install the script `simutator`.

## Tests

Run `tox` to run all the tests.


## Reproducibility

Both commands below are by default non-deterministic.
However, both have a `--seed` option, which you should use if you want to
make reproducible runs from the same input data. For example, run everything
with `--seed 42`.

## Make mutated genomes

The general command is
```
simutator mutate_fasta [options] in.fasta out
```

At least one of these options must be used: `--snps`, `--dels`, `--ins`,
`--complex`.

### Example 1 - simple

Simulate a SNP every 100bp with

```
simutator mutate_fasta --snps 100 in.fasta out
```

This will write three output files:
1. `out.snp.dist-100.fa` - a mutated version of `in.fasta`, with SNPs added
2. `out.snp.dist-100.original.vcf` - VCF file of variants between the original
   and mutated FASTA files, where the VCF reference is the _original_ FASTA.
3. `out.snp.dist-100.mutated.vcf` - VCF file of variants between the original
   and mutated FASTA files, where the VCF reference is the _mutated_ FASTA.


### Example 2 - reasonably simple

Independently simulate a SNP every 100bp, a SNP every 200bp, and a
3bp deletion every 2kb with

```
simutator mutate_fasta --snps 100,200 --dels 2000:3 in.fasta out
```

That will write a new FASTA and two VCF files for the SNPs every 100bp, and
similarly for SNPs every 200bp and for the deletions. ie each of the three
simulations are completely independent of each other and result in
three new files each (a FASTA and two VCF files).

### Example 3 - expert

Make:
1. SNP every 100bp
2. SNP every 2000bp
3. deletion of length 3 every 2kb
4. deletion of length 4 every 42kb
5. insertion of length 7 every 10kb
6. "complex variant" every 2kb where a region of length 20bp gets
  3 snps, 1 insertion and 2 deletions up to length 4
7. "complex variant" every 10kb where a region of length 50bp gets
  10 snps, 2 insertions and 3 deletions up to length 5

```
simutator mutate_fasta \
  --snps 100,2000
  --dels 200:3,42000:4 \
  --ins 10000:7 \
  --complex 2000:20:3:1:2:4,10000:50:10:2:3:5 \
  in.fasta out
```

That will make 7 sets of output files, one set of three files for each of
the mutation types in the above list. A FASTA of the mutated genome, a VCF
where the reference is the original genome, and a VCF where the reference is
the mutated genome.

## Make simulated reads

Requires `art_illumina` to be in your `$PATH` (see the dependencies section above).
If you use this, please cite "ART: a next-generation sequencing read simulator",
Huang et al 2011.

The command is

```
simutator simulate_reads in.fasta out [options]
```

This is simply a wrapper around `art_illumina` that can generate combinations
of sequencing machines, read lengths, read depths, and fragment lengths.

For example, to generate all combinations of read lengths 100 and 150,
and read depths 1X and 2X:

```
simutator in.fasta out --read_length 100 150 --read_depth 1 2
```

That made four sets of read pairs, which means eight (gzipped FASTQ) files:
```
out.HS25.100.1.500.25.1.fq.gz
out.HS25.100.1.500.25.2.fq.gz
out.HS25.100.2.500.25.1.fq.gz
out.HS25.100.2.500.25.2.fq.gz
out.HS25.150.1.500.25.1.fq.gz
out.HS25.150.1.500.25.2.fq.gz
out.HS25.150.2.500.25.1.fq.gz
out.HS25.150.2.500.25.2.fq.gz
```
and a JSON file `out.json` with the details of each set of reads:
```json
[
  {
    "fastq1": "out.HS25.100.1.500.25.1.fq.gz",
    "fastq2": "out.HS25.100.1.500.25.2.fq.gz",
    "fragment_length": 500,
    "fragment_length_sd": 25,
    "machine": "HS25",
    "read_depth": 1,
    "read_length": 100
  },
  {
    "fastq1": "out.HS25.100.2.500.25.1.fq.gz",
    "fastq2": "out.HS25.100.2.500.25.2.fq.gz",
    "fragment_length": 500,
    "fragment_length_sd": 25,
    "machine": "HS25",
    "read_depth": 2,
    "read_length": 100
  },
  {
    "fastq1": "out.HS25.150.1.500.25.1.fq.gz",
    "fastq2": "out.HS25.150.1.500.25.2.fq.gz",
    "fragment_length": 500,
    "fragment_length_sd": 25,
    "machine": "HS25",
    "read_depth": 1,
    "read_length": 150
  },
  {
    "fastq1": "out.HS25.150.2.500.25.1.fq.gz",
    "fastq2": "out.HS25.150.2.500.25.2.fq.gz",
    "fragment_length": 500,
    "fragment_length_sd": 25,
    "machine": "HS25",
    "read_depth": 2,
    "read_length": 150
  }
]
```
