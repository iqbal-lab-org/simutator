import json

from simutator import simulate_reads


def run(options):
    data = simulate_reads.iterative_simulate_reads(
        options.fasta_in,
        options.outprefix,
        options.machine,
        options.read_length,
        options.read_depth,
        options.fragment_length,
        options.fragment_length_sd,
        random_seed=options.seed,
    )
    with open(options.outprefix + ".json", "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)
