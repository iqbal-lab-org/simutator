import logging

from simutator import genome_mutator


def _parse_indels_option_string(s):
    distances_and_lengths = []
    for x in s.split(","):
        dist, length = x.split(":")
        distances_and_lengths.append({"dist": int(dist), "len": int(length)})
    return distances_and_lengths


def _parse_complex_option_string(s):
    complex_vars = []
    for x in s.split(","):
        dist, length, snps, ins, dels, max_indel = [int(y) for y in x.split(":")]
        complex_vars.append(
            {
                "dist": dist,
                "len": length,
                "snp": snps,
                "ins": ins,
                "del": dels,
                "max_indel_len": max_indel,
            }
        )
    return complex_vars


def mutations_from_options(options):
    mutations = {}
    if options.snps is not None:
        try:
            snp_distances = [{"dist": int(x)} for x in options.snps.split(",")]
        except:
            raise ValueError(f"Cannot parse --snps option: '{options.snps}'")
        mutations["snp"] = snp_distances

    if options.ins is not None:
        try:
            mutations["insertion"] = _parse_indels_option_string(options.ins)
        except:
            raise ValueError(f"Cannot parse --ins option: '{options.ins}'")

    if options.dels is not None:
        try:
            mutations["deletion"] = _parse_indels_option_string(options.dels)
        except:
            raise ValueError(f"Cannot parse --dels option: '{options.dels}'")

    if options.complex is not None:
        try:
            mutations["complex"] = _parse_complex_option_string(options.complex)
        except:
            raise ValueError(f"Cannot parse --complex option: '{options.complex}'")

    if len(mutations) == 0:
        raise RuntimeError(
            "Must use at least one of the options --snps, --dels, --ins, --complex"
        )

    return mutations


def run_all_mutations(fasta_in, outprefix, mutations, seed=None):
    for mutation_type, mutations_list in mutations.items():
        for mutation in mutations_list:
            logging.info(
                f"Simulating mutations of type '{mutation_type}' with parameters {mutation}"
            )
            if mutation_type == "snp":
                mutator = genome_mutator.SnpMutator(mutation["dist"], seed=seed)
            elif mutation_type == "insertion":
                mutator = genome_mutator.InsertionMutator(
                    mutation["dist"], mutation["len"], seed=seed
                )
            elif mutation_type == "deletion":
                mutator = genome_mutator.DeletionMutator(
                    mutation["dist"], mutation["len"], seed=seed
                )
            elif mutation_type == "complex":
                mutator = genome_mutator.ComplexMutator(
                    mutation["dist"],
                    mutation["len"],
                    mutation["snp"],
                    mutation["del"],
                    mutation["ins"],
                    mutation["max_indel_len"],
                    seed=seed,
                )

            this_prefix = f"{outprefix}.{mutation_type}." + ".".join(
                [k + "-" + str(v) for k, v in sorted(mutation.items())]
            )
            mutator.mutate_fasta_file(
                fasta_in,
                f"{this_prefix}.fa",
                f"{this_prefix}.original.vcf",
                f"{this_prefix}.mutated.vcf",
            )
