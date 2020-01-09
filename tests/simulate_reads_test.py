import os
import pytest

import pyfastaq

from simutator import simulate_reads, utils


def test_simulate_illumina_paired_reads_from_fasta():
    """test simulate_illumina_paired_reads_from_fasta"""
    #  Just test if we can run ART. Make a fake genome and then simulate
    # some reads from it
    tmp_outprefix = "tmp.simulate_illumina_paired_reads_from_fasta"
    utils.syscall(f"rm -rf {tmp_outprefix}.*")
    p = utils.syscall("ls -lh")
    print(p.stdout)
    tmp_ref = f"{tmp_outprefix}.ref.fa"
    pyfastaq.tasks.make_random_contigs(2, 2000, tmp_ref)
    outfiles = simulate_reads.simulate_illumina_paired_reads_from_fasta(
        tmp_ref, tmp_outprefix
    )
    os.unlink(tmp_ref)
    for filename in outfiles:
        assert os.path.exists(filename)
        os.unlink(filename)
