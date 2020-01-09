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
