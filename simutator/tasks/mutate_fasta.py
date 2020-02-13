from simutator import batch_genome_mutator


def run(options):
    mutations = batch_genome_mutator.mutations_from_options(options)
    batch_genome_mutator.run_all_mutations(
        options.fasta_in, options.outprefix, mutations, seed=options.seed
    )
