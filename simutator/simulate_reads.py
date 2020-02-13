import itertools
import logging
import os
import shutil
import tempfile

from simutator import utils

# This uses ART to simulated reads. Get it like this:
#   wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier20160605linux64tgz.tgz
#   tar xf artbinmountrainier20160605linux64tgz.tgz
# The executable is:
#   $PWD/art_bin_MountRainier/art_illumina
def simulate_illumina_paired_reads_from_fasta(
    ref_fasta,
    outprefix,
    sequencing_machine="HS25",
    read_length=150,
    read_depth=50,
    mean_fragment_length=500,
    fragment_length_sd=25,
    random_seed=42,
):
    """Simulates Illumina paired end reads using ART.
    Returns tuple (forward reads filename, reverse reads filename)"""
    if shutil.which("art_illumina") is None:
        raise RuntimeError("art_illumina not found in PATH. Cannot continue")

    seed_string = "" if random_seed is None else "--rndSeed " + str(random_seed)
    tmpdir = tempfile.mkdtemp(prefix=outprefix + ".", dir=os.getcwd())
    tmp_prefix = os.path.join(tmpdir, "out")

    command = " ".join(
        [
            "art_illumina",
            "--in",
            ref_fasta,
            "--out",
            tmp_prefix,
            "--noALN",  # do not output alignment file
            "--seqSys",
            sequencing_machine,
            "--len",
            str(read_length),
            "--fcov",
            str(read_depth),
            "--mflen",
            str(mean_fragment_length),
            "--sdev",
            str(fragment_length_sd),
            seed_string,
        ]
    ).rstrip()

    utils.syscall(command)
    reads_files = []

    for i in ("1", "2"):
        final_reads_file = outprefix + "." + i + ".fq"
        os.rename(tmp_prefix + i + ".fq", final_reads_file)
        utils.syscall("gzip -9 " + final_reads_file)
        reads_files.append(f"{final_reads_file}.gz")

    os.rmdir(tmpdir)
    return tuple(reads_files)


def iterative_simulate_reads(
    ref_fasta,
    outprefix,
    sequencing_machines,
    read_lengths,
    read_depths,
    fragment_lengths,
    fragment_length_sd,
    random_seed=42,
):
    files = []

    for machine, read_len, depth, frag_len in itertools.product(
        sequencing_machines, read_lengths, read_depths, fragment_lengths
    ):
        this_prefix = (
            f"{outprefix}.{machine}.{read_len}.{depth}.{frag_len}.{fragment_length_sd}"
        )
        logging.info(
            "Simulate reads. ref={ref_fasta}, machine={machine}, read length={read_len}, read depth={depth}, fragment length={frag_len}, fragment length sd={fragment_length_sd}"
        )

        reads_files = simulate_illumina_paired_reads_from_fasta(
            ref_fasta,
            this_prefix,
            sequencing_machine=machine,
            read_length=read_len,
            read_depth=depth,
            mean_fragment_length=frag_len,
            fragment_length_sd=fragment_length_sd,
            random_seed=random_seed,
        )

        files.append(
            {
                "fastq1": reads_files[0],
                "fastq2": reads_files[1],
                "machine": machine,
                "read_length": read_len,
                "read_depth": depth,
                "fragment_length": frag_len,
                "fragment_length_sd": fragment_length_sd,
            }
        )

    return files
