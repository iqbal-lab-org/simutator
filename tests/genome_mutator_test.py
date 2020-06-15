import filecmp
import os
import pytest

import pyfastaq

from simutator import genome_mutator

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "genome_mutator")


def test_SnpMutator_mutate_sequence():
    # Set the seed so that the random mutations will be the same every time
    # the test runs
    mutator = genome_mutator.SnpMutator(3, seed=42)
    original_seq = "AGTAGGCAG"
    sequence = pyfastaq.sequences.Fasta("name", original_seq)
    got_mutations, got_sequence = mutator.mutate_sequence(sequence)
    assert sequence.seq == original_seq
    assert got_sequence == "AGGAGACAG"
    expect_mutations = [
        genome_mutator.Mutation(2, 2, "T", "G"),
        genome_mutator.Mutation(5, 5, "G", "A"),
    ]
    assert got_mutations == expect_mutations


def test_SnpMutator_mutate_fasta_file():
    infile = os.path.join(data_dir, "SnpMutator_mutate_fasta.in.fa")
    expected_fa = os.path.join(data_dir, "SnpMutator_mutate_fasta.out.fa")
    expected_vcf_ref = os.path.join(data_dir, "SnpMutator_mutate_fasta.out.ref.vcf")
    expected_vcf_mutated = os.path.join(
        data_dir, "SnpMutator_mutate_fasta.out.mutated.vcf"
    )
    mutator = genome_mutator.SnpMutator(30, seed=42)
    tmp_out_fa = "tmp.SnpMutator_mutate_fasta.out.fa"
    tmp_out_vcf_ref = "tmp.SnpMutator_mutate_fasta.out.ref.vcf"
    tmp_out_vcf_mutated = "tmp.SnpMutator_mutate_fasta.out.mutated.vcf"
    mutator.mutate_fasta_file(infile, tmp_out_fa, tmp_out_vcf_ref, tmp_out_vcf_mutated)
    assert filecmp.cmp(tmp_out_fa, expected_fa, shallow=False)
    assert filecmp.cmp(tmp_out_vcf_ref, expected_vcf_ref, shallow=False)
    assert filecmp.cmp(tmp_out_vcf_mutated, expected_vcf_mutated, shallow=False)
    os.unlink(tmp_out_fa)
    os.unlink(tmp_out_vcf_ref)
    os.unlink(tmp_out_vcf_mutated)


def test_DeletionMutator_mutate_sequence_del_length_1():
    mutator = genome_mutator.DeletionMutator(3, 1)
    original_seq = "1234567890ABCDE"
    sequence = pyfastaq.sequences.Fasta("name", original_seq)
    got_mutations, got_sequence = mutator.mutate_sequence(sequence)
    assert sequence.seq == original_seq
    assert got_sequence == "1245780ACDE"
    expect_mutations = [
        genome_mutator.Mutation(1, 1, "23", "2"),
        genome_mutator.Mutation(4, 3, "56", "5"),
        genome_mutator.Mutation(7, 5, "89", "8"),
        genome_mutator.Mutation(10, 7, "AB", "A"),
    ]
    assert got_mutations == expect_mutations


def test_DeletionMutator_mutate_sequence_del_length_2():
    mutator = genome_mutator.DeletionMutator(5, 2)
    original_seq = "1234567890ABCDEF"
    sequence = pyfastaq.sequences.Fasta("name", original_seq)
    got_mutations, got_sequence = mutator.mutate_sequence(sequence)
    assert sequence.seq == original_seq
    assert got_sequence == "12347890CDEF"
    expect_mutations = [
        genome_mutator.Mutation(3, 3, "456", "4"),
        genome_mutator.Mutation(9, 7, "0AB", "0"),
    ]
    assert got_mutations == expect_mutations

def test_DeletionMutator_mutate_fasta_file():
    infile = os.path.join(data_dir, "DeletionMutator_mutate_fasta.in.fa")
    expected_fa = os.path.join(data_dir, "DeletionMutator_mutate_fasta.out.fa")
    expected_vcf_ref = os.path.join(
        data_dir, "DeletionMutator_mutate_fasta.out.ref.vcf"
    )
    expected_vcf_mutated = os.path.join(
        data_dir, "DeletionMutator_mutate_fasta.out.mutated.vcf"
    )
    mutator = genome_mutator.DeletionMutator(30, 1, seed=42)
    tmp_out_fa = "tmp.DeletionMutator_mutate_fasta.out.fa"
    tmp_out_vcf_ref = "tmp.DeletionMutator_mutate_fasta.out.ref.vcf"
    tmp_out_vcf_mutated = "tmp.DeletionMutator_mutate_fasta.out.mutated.vcf"
    mutator.mutate_fasta_file(infile, tmp_out_fa, tmp_out_vcf_ref, tmp_out_vcf_mutated)
    assert filecmp.cmp(tmp_out_fa, expected_fa, shallow=False)
    assert filecmp.cmp(tmp_out_vcf_ref, expected_vcf_ref, shallow=False)
    assert filecmp.cmp(tmp_out_vcf_mutated, expected_vcf_mutated, shallow=False)
    os.unlink(tmp_out_fa)
    os.unlink(tmp_out_vcf_ref)
    os.unlink(tmp_out_vcf_mutated)


def test_InsertionMutator_mutate_sequence():
    mutator = genome_mutator.InsertionMutator(4, 1, seed=43)
    original_seq = "0123456789XYZ"
    sequence = pyfastaq.sequences.Fasta("name", original_seq)
    got_mutations, got_sequence = mutator.mutate_sequence(sequence)
    assert sequence.seq == original_seq
    assert got_sequence == "0123A4567G89XYZ"
    expect_mutations = [
        genome_mutator.Mutation(3, 3, "3", "3A"),
        genome_mutator.Mutation(7, 8, "7", "7G"),
    ]
    assert got_mutations == expect_mutations


def test_InsertionMutator_mutate_fasta_file():
    infile = os.path.join(data_dir, "InsertionMutator_mutate_fasta.in.fa")
    expected_fa = os.path.join(data_dir, "InsertionMutator_mutate_fasta.out.fa")
    expected_vcf_ref = os.path.join(
        data_dir, "InsertionMutator_mutate_fasta.out.ref.vcf"
    )
    expected_vcf_mutated = os.path.join(
        data_dir, "InsertionMutator_mutate_fasta.out.mutated.vcf"
    )
    mutator = genome_mutator.InsertionMutator(30, 1, seed=42)
    tmp_out_fa = "tmp.InsertionMutator_mutate_fasta.out.fa"
    tmp_out_vcf_ref = "tmp.InsertionMutator_mutate_fasta.out.ref.vcf"
    tmp_out_vcf_mutated = "tmp.InsertionMutator_mutate_fasta.out.mutated.vcf"
    mutator.mutate_fasta_file(infile, tmp_out_fa, tmp_out_vcf_ref, tmp_out_vcf_mutated)
    assert filecmp.cmp(tmp_out_fa, expected_fa, shallow=False)
    assert filecmp.cmp(tmp_out_vcf_ref, expected_vcf_ref, shallow=False)
    assert filecmp.cmp(tmp_out_vcf_mutated, expected_vcf_mutated, shallow=False)
    os.unlink(tmp_out_fa)
    os.unlink(tmp_out_vcf_ref)
    os.unlink(tmp_out_vcf_mutated)


def test_ComplexMutator_mutate_sequence():
    mutator = genome_mutator.ComplexMutator(20, 10, 2, 2, 1, 1, seed=42)
    original_seq = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgh"
    sequence = pyfastaq.sequences.Fasta("name", original_seq)
    got_mutations, got_sequence = mutator.mutate_sequence(sequence)
    print(got_sequence)
    assert sequence.seq == original_seq
    assert got_sequence == "abcdefghijklmnopqrstAuwAyabTdefghijklmnTqCtTuvwxyzabcdefgh"
    expect_mutations = [
        genome_mutator.Mutation(19, 19, "tuvwxyzabc", "tAuwAyabT"),
        genome_mutator.Mutation(39, 38, "nopqrstuvw", "nTqCtTuvw"),
    ]
    assert got_mutations == expect_mutations


def test_ComplexMutator_mutate_fasta_file():
    infile = os.path.join(data_dir, "ComplexMutator_mutate_fasta.in.fa")
    expected_fa = os.path.join(data_dir, "ComplexMutator_mutate_fasta.out.fa")
    expected_vcf_ref = os.path.join(data_dir, "ComplexMutator_mutate_fasta.out.ref.vcf")
    expected_vcf_mutated = os.path.join(
        data_dir, "ComplexMutator_mutate_fasta.out.mutated.vcf"
    )
    mutator = genome_mutator.ComplexMutator(30, 10, 2, 2, 1, 2, seed=42)
    tmp_out_fa = "tmp.ComplexMutator_mutate_fasta.out.fa"
    tmp_out_vcf_ref = "tmp.ComplexMutator_mutate_fasta.out.ref.vcf"
    tmp_out_vcf_mutated = "tmp.ComplexMutator_mutate_fasta.out.mutated.vcf"
    mutator.mutate_fasta_file(infile, tmp_out_fa, tmp_out_vcf_ref, tmp_out_vcf_mutated)
    assert filecmp.cmp(tmp_out_fa, expected_fa, shallow=False)
    assert filecmp.cmp(tmp_out_vcf_ref, expected_vcf_ref, shallow=False)
    assert filecmp.cmp(tmp_out_vcf_mutated, expected_vcf_mutated, shallow=False)
    os.unlink(tmp_out_fa)
    os.unlink(tmp_out_vcf_ref)
    os.unlink(tmp_out_vcf_mutated)
