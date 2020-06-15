import abc
import collections
from random import Random

import pyfastaq

global random  # to keep seeding consistent
random = Random()

Mutation = collections.namedtuple(
    "Mutation", ["original_position", "new_position", "original_seq", "new_seq"]
)
acgt = {"A", "C", "G", "T"}


class GenomeMutator(metaclass=abc.ABCMeta):
    def __init__(self, distance_between_mutations, seed=None):
        self.distance_between_mutations = distance_between_mutations
        if seed is not None:
            global random
            random = Random(seed)

    def _vcf_source_prefix(self, mutated_genome=False):
        if mutated_genome:
            return "##source=simutator, ref in this file is mutated genome. Mutations added:"
        else:
            return "##source=simutator, ref in this file is original genome. Mutations added:"

    @abc.abstractmethod
    def _mutation_description_string(self):
        pass

    def _vcf_source_line(self, mutated_genome=False):
        return (
            self._vcf_source_prefix(mutated_genome=mutated_genome)
            + " "
            + self._mutation_description_string()
        )

    def _write_vcf_header(self, filehandle, seq_lengths, mutated_genome=False):
        print("##fileformat=VCFv4.2", file=filehandle)
        print(self._vcf_source_line(mutated_genome=mutated_genome), file=filehandle)
        for name, length in sorted(seq_lengths.items()):
            print(f"##contig=<ID={name},length={length}>", file=filehandle)
        print(
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "sample",
            sep="\t",
            file=filehandle,
        )

    @abc.abstractmethod
    def mutate_sequence(self, sequence):
        pass

    def mutate_fasta_file(
        self, fasta_in, fasta_out, vcf_out_wrt_original_seq, vcf_out_wrt_mutated_seq
    ):
        file_reader = pyfastaq.sequences.file_reader(fasta_in)
        original_seq_lengths = {}
        mutated_seq_lengths = {}
        all_mutations = {}

        with open(fasta_out, "w") as f_fasta:
            for sequence in file_reader:
                mutations, mutated_seq = self.mutate_sequence(sequence)
                mutated_seq = pyfastaq.sequences.Fasta(
                    sequence.id + "__simutator__" + self._mutation_description_string(),
                    mutated_seq,
                )
                print(mutated_seq, file=f_fasta)
                original_seq_lengths[sequence.id] = len(sequence)
                mutated_seq_lengths[mutated_seq.id] = len(mutated_seq)
                all_mutations[(sequence.id, mutated_seq.id)] = mutations

        with open(vcf_out_wrt_original_seq, "w") as f_vcf_original, open(
            vcf_out_wrt_mutated_seq, "w"
        ) as f_vcf_mutated:
            self._write_vcf_header(
                f_vcf_original, original_seq_lengths, mutated_genome=False
            )
            self._write_vcf_header(
                f_vcf_mutated, mutated_seq_lengths, mutated_genome=True
            )

            for seq_id, mutated_seq_id in sorted(all_mutations):
                for mutation in all_mutations[(seq_id, mutated_seq_id)]:
                    print(
                        mutated_seq_id,
                        mutation.new_position + 1,
                        ".",
                        mutation.new_seq,
                        mutation.original_seq,
                        ".",
                        "PASS",
                        ".",
                        "GT",
                        "1/1",
                        sep="\t",
                        file=f_vcf_mutated,
                    )
                    print(
                        seq_id,
                        mutation.original_position + 1,
                        ".",
                        mutation.original_seq,
                        mutation.new_seq,
                        ".",
                        "PASS",
                        ".",
                        "GT",
                        "1/1",
                        sep="\t",
                        file=f_vcf_original,
                    )

    def _get_snp_variant(self, ref_nucleotide):
        global random
        return random.choice(sorted(list(acgt.difference({ref_nucleotide}))))


class SnpMutator(GenomeMutator):
    def __init__(self, distance_between_snps, seed=None):
        super().__init__(distance_between_snps, seed=seed)

    def _mutation_description_string(self):
        return f"SNP_every_{self.distance_between_mutations}"

    def mutate_sequence(self, sequence):
        mutations = []
        new_sequence = list(sequence)
        for i in range(
            self.distance_between_mutations - 1,
            len(sequence) - self.distance_between_mutations,
            self.distance_between_mutations,
        ):
            old_nucleotide = new_sequence[i].upper()
            new_sequence[i] = self._get_snp_variant(old_nucleotide)
            mutations.append(Mutation(i, i, old_nucleotide, new_sequence[i]))

        mutated_seq = "".join(new_sequence)
        return mutations, mutated_seq


class DeletionMutator(GenomeMutator):
    def __init__(self, distance_between_deletions, deletion_length, seed=None):
        super().__init__(distance_between_deletions, seed=seed)
        self.deletion_length = deletion_length

    def _mutation_description_string(self):
        return (
            f"DEL_length_{self.deletion_length}_every_{self.distance_between_mutations}"
        )

    def mutate_sequence(self, sequence):
        mutations = []
        current_position = self.distance_between_mutations - 1
        deleted_nucleotides = 0
        new_sequence = [sequence[:current_position]]

        while current_position < len(sequence) - self.distance_between_mutations:
            deletion_end_position = current_position + self.deletion_length - 1
            next_start_position = (
                deletion_end_position + self.distance_between_mutations
            )
            new_sequence.append(
                sequence.seq[deletion_end_position + 1 : next_start_position]
            )
            mutations.append(
                Mutation(
                    current_position - 1,
                    current_position - deleted_nucleotides - 1,
                    sequence.seq[current_position - 1 : deletion_end_position + 1],
                    sequence.seq[current_position - 1],
                )
            )
            current_position = next_start_position
            deleted_nucleotides += self.deletion_length

        new_sequence.append(sequence.seq[current_position:])
        mutated_seq = "".join(new_sequence)
        return mutations, mutated_seq


class InsertionMutator(GenomeMutator):
    def __init__(self, distance_between_insertions, insertion_length, seed=None):
        super().__init__(distance_between_insertions, seed=seed)
        self.insertion_length = insertion_length

    def _mutation_description_string(self):
        return f"INS_length_{self.insertion_length}_every_{self.distance_between_mutations}"

    def mutate_sequence(self, sequence):
        mutations = []
        current_position = self.distance_between_mutations
        inserted_nucleotides = 0
        new_sequence = [sequence.seq[:current_position]]

        while current_position < len(sequence) - self.distance_between_mutations:
            insertion_seq = "".join(
                [
                    random.choice(["A", "C", "G", "T"])
                    for _ in range(self.insertion_length)
                ]
            )
            new_sequence.append(insertion_seq)
            new_sequence.append(
                sequence.seq[
                    current_position : current_position
                    + self.distance_between_mutations
                ]
            )
            mutations.append(
                Mutation(
                    current_position - 1,
                    current_position + inserted_nucleotides - 1,
                    sequence.seq[current_position - 1],
                    sequence.seq[current_position - 1] + insertion_seq,
                )
            )
            current_position += self.distance_between_mutations
            inserted_nucleotides += len(insertion_seq)

        new_sequence.append(sequence.seq[current_position:])
        mutated_seq = "".join(new_sequence)
        return mutations, mutated_seq


class ComplexMutator(GenomeMutator):
    def __init__(
        self,
        distance_between_clusters,
        cluster_length,
        snps_per_cluster,
        dels_per_cluster,
        ins_per_cluster,
        max_indel_length,
        seed=None,
    ):
        super().__init__(distance_between_clusters, seed=seed)
        self.cluster_length = cluster_length
        self.snps_per_cluster = snps_per_cluster
        self.dels_per_cluster = dels_per_cluster
        self.ins_per_cluster = ins_per_cluster
        self.max_indel_length = max_indel_length

    def _mutation_description_string(self):
        return "_".join(
            [
                f"COMPLEX_length_{self.cluster_length}",
                f"every_{self.distance_between_mutations}",
                f"snps_{self.snps_per_cluster}",
                f"del_{self.dels_per_cluster}",
                f"ins_{self.ins_per_cluster}",
                f"maxindel_{self.max_indel_length}",
            ]
        )

    def _add_cluster_of_variants_to_sequence(
        self, sequence, deletion_lengths, insertion_lengths
    ):
        total_variations = (
            self.snps_per_cluster + len(deletion_lengths) + len(insertion_lengths)
        )
        variant_positions = random.sample(range(1, len(sequence)), total_variations)
        snp_positions = sorted(variant_positions[: self.snps_per_cluster])
        deletion_positions = variant_positions[
            self.snps_per_cluster : self.snps_per_cluster + len(deletion_lengths)
        ]
        insertion_positions = variant_positions[
            self.snps_per_cluster + len(deletion_lengths) :
        ]
        indels = {
            pos: (length, "ins")
            for (pos, length) in zip(insertion_positions, insertion_lengths)
        }
        indels.update(
            {
                pos: (length, "del")
                for (pos, length) in zip(deletion_positions, deletion_lengths)
            }
        )
        nucleotides_list = list(sequence)

        for snp_position in snp_positions:
            nucleotides_list[snp_position] = self._get_snp_variant(
                nucleotides_list[snp_position]
            )

        position_offset = 0

        for indel_position in sorted(indels):
            offset_position = indel_position + position_offset
            if offset_position > len(nucleotides_list) - 1:
                continue

            indel_length, ins_or_del = indels[indel_position]
            if ins_or_del == "ins":
                nucleotides_list[offset_position:offset_position] = [
                    random.choice(["A", "C", "G", "T"]) for _ in range(indel_length)
                ]
                position_offset += indel_length
            else:
                assert ins_or_del == "del"
                del nucleotides_list[offset_position : offset_position + indel_length]
                position_offset -= indel_length

        return "".join(nucleotides_list)

    def mutate_sequence(self, sequence):
        new_sequence = [sequence.seq[: self.distance_between_mutations - 1]]
        mutations = []
        cluster_start = None
        inserted_bases = 0

        for cluster_start in range(
            self.distance_between_mutations - 1,
            len(sequence) - self.distance_between_mutations,
            self.distance_between_mutations,
        ):
            deletion_lengths = [
                random.randint(1, self.max_indel_length)
                for _ in range(self.dels_per_cluster)
            ]
            insertion_lengths = [
                random.randint(1, self.max_indel_length)
                for _ in range(self.ins_per_cluster)
            ]
            original_cluster_seq = sequence.seq[
                cluster_start : cluster_start + self.cluster_length
            ]
            variant_seq = self._add_cluster_of_variants_to_sequence(
                original_cluster_seq, deletion_lengths, insertion_lengths
            )
            new_sequence.append(
                variant_seq
                + sequence.seq[
                    cluster_start
                    + self.cluster_length : cluster_start
                    + self.distance_between_mutations
                ]
            )
            mutations.append(
                Mutation(
                    cluster_start,
                    cluster_start + inserted_bases,
                    original_cluster_seq,
                    variant_seq,
                )
            )
            inserted_bases += len(variant_seq) - len(original_cluster_seq)

        if cluster_start is not None:
            cluster_start += self.distance_between_mutations
            new_sequence.append(sequence.seq[cluster_start:])

        mutated_seq = "".join(new_sequence)
        return mutations, mutated_seq
