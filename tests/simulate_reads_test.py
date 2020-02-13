import os
import pytest
import shutil

import pyfastaq

from simutator import simulate_reads, utils


def test_simulate_illumina_paired_reads_from_fasta():
    """test simulate_illumina_paired_reads_from_fasta"""
    #  Just test if we can run ART. Make a fake genome and then simulate
    # some reads from it
    tmp_outprefix = "tmp.simulate_illumina_paired_reads_from_fasta"
    utils.syscall(f"rm -rf {tmp_outprefix}.*")
    tmp_ref = f"{tmp_outprefix}.ref.fa"
    pyfastaq.tasks.make_random_contigs(2, 2000, tmp_ref)
    outfiles = simulate_reads.simulate_illumina_paired_reads_from_fasta(
        tmp_ref, tmp_outprefix
    )
    os.unlink(tmp_ref)
    for filename in outfiles:
        assert os.path.exists(filename)
        os.unlink(filename)


def test_iterative_simulate_reads():
    tmpdir = "tmp.iterative_simulate_reads"
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.mkdir(tmpdir)
    tmp_ref = os.path.join(tmpdir, "ref.fa")
    outprefix = os.path.join(tmpdir, "out")
    pyfastaq.tasks.make_random_contigs(1, 2000, tmp_ref)
    got = simulate_reads.iterative_simulate_reads(
        tmp_ref, outprefix, ["HS25"], [50, 100], [1, 2], [300], 10, random_seed=42
    )
    expect = [
        {
            "fastq1": "tmp.iterative_simulate_reads/out.HS25.50.1.300.10.1.fq.gz",
            "fastq2": "tmp.iterative_simulate_reads/out.HS25.50.1.300.10.2.fq.gz",
            "fragment_length": 300,
            "fragment_length_sd": 10,
            "machine": "HS25",
            "read_depth": 1,
            "read_length": 50,
        },
        {
            "fastq1": "tmp.iterative_simulate_reads/out.HS25.50.2.300.10.1.fq.gz",
            "fastq2": "tmp.iterative_simulate_reads/out.HS25.50.2.300.10.2.fq.gz",
            "fragment_length": 300,
            "fragment_length_sd": 10,
            "machine": "HS25",
            "read_depth": 2,
            "read_length": 50,
        },
        {
            "fastq1": "tmp.iterative_simulate_reads/out.HS25.100.1.300.10.1.fq.gz",
            "fastq2": "tmp.iterative_simulate_reads/out.HS25.100.1.300.10.2.fq.gz",
            "fragment_length": 300,
            "fragment_length_sd": 10,
            "machine": "HS25",
            "read_depth": 1,
            "read_length": 100,
        },
        {
            "fastq1": "tmp.iterative_simulate_reads/out.HS25.100.2.300.10.1.fq.gz",
            "fastq2": "tmp.iterative_simulate_reads/out.HS25.100.2.300.10.2.fq.gz",
            "fragment_length": 300,
            "fragment_length_sd": 10,
            "machine": "HS25",
            "read_depth": 2,
            "read_length": 100,
        },
    ]
    assert got == expect
    for d in expect:
        assert os.path.exists(d["fastq1"])
        assert os.path.exists(d["fastq2"])
    shutil.rmtree(tmpdir)
