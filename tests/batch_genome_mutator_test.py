import os
import pytest
import shutil

from simutator import batch_genome_mutator

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "batch_genome_mutator")


def test_parse_indels_option_string():
    with pytest.raises(ValueError):
        batch_genome_mutator._parse_indels_option_string("this_is_not_even_close")
    with pytest.raises(ValueError):
        batch_genome_mutator._parse_indels_option_string("1:not_a_number")

    expect = [{"dist": 100, "len": 3}, {"dist": 500, "len": 10}]
    got = batch_genome_mutator._parse_indels_option_string("100:3,500:10")
    assert got == expect


def test_parse_complex_option_string():
    with pytest.raises(ValueError):
        batch_genome_mutator._parse_complex_option_string("totally_unexpected")

    got = batch_genome_mutator._parse_complex_option_string(
        "1000:10:1:2:3:4,500:20:2:3:4:5"
    )
    expect = [
        {"dist": 1000, "len": 10, "snp": 1, "ins": 2, "del": 3, "max_indel_len": 4},
        {"dist": 500, "len": 20, "snp": 2, "ins": 3, "del": 4, "max_indel_len": 5},
    ]
    assert got == expect


def test_run_all_mutations():
    infile = os.path.join(data_dir, "run_all_mutations.fa")
    mutations = {
        "snp": [{"dist": 200}, {"dist": 300}],
        "ins": [{"dist": 200, "len": 10}],
        "del": [{"dist": 250, "len": 5}],
        "complex": [
            {"dist": 500, "len": 20, "snp": 2, "ins": 3, "del": 4, "max_indel_len": 5}
        ],
    }
    outdir = "tmp.run_all_mutations"
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.mkdir(outdir)
    outprefix = os.path.join(outdir, "out")
    batch_genome_mutator.run_all_mutations(infile, outprefix, mutations)

    expect_prefixes = [
        outprefix + "." + x
        for x in [
            "complex.del-4.dist-500.ins-3.len-20.max_indel_len-5.snp-2",
            "del.dist-250.len-5",
            "ins.dist-200.len-10",
            "snp.dist-200",
            "snp.dist-300",
        ]
    ]

    for prefix in expect_prefixes:
        for suffix in "fa", "mutated.vcf", "original.vcf":
            assert os.path.exists(f"{prefix}.{suffix}")

    shutil.rmtree(outdir)
