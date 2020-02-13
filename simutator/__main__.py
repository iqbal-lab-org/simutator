#!/usr/bin/env python3

import argparse
import logging
import simutator

def main(args=None):
    parser = argparse.ArgumentParser(
        prog="simutator",
        usage="simutator <command> <options>",
        description="simutator: simulate mutations in genomes"
    )

    parser.add_argument("--version", action="version", version=simutator.__version__)
    parser.add_argument("--debug", help="Debug mode", action="store_true")
    subparsers = parser.add_subparsers(title="Available commands", help="", metavar="")

    # ----------------------- mutate_fasta ----------------------------------------
    subparser_mutate_fasta = subparsers.add_parser(
        "mutate_fasta",
        help="Mutate a FASTA file",
        usage="simutator mutate_fasta [options] <in.fasta> <outprefix>",
        description="Mutate a FASTA file",
    )

    subparser_mutate_fasta.add_argument(
        "--seed",
        help="Seed for random number generator. Use this option for reproducibility, otherwise Python's default seeding is used",
        type=int,
        metavar="INT",
    )

    subparser_mutate_fasta.add_argument(
        "--snps",
        help="Comma-separated list of distances between SNPs",
        metavar="INT[,INT,...]"
    )

    subparser_mutate_fasta.add_argument(
        "--dels",
        help="Comma-separated list of <distance between deletions>:<deletion lengths>",
        metavar="INT1:INT2[,INT3:INT4,...]"
    )

    subparser_mutate_fasta.add_argument(
        "--ins",
        help="Comma-separated list of <distance between insertions>:<insertion lengths>",
        metavar="INT1:INT2[,INT3:INT4,...]"
    )

    subparser_mutate_fasta.add_argument(
        "--complex",
        help="Comma-separated list of dist:len:s:ins:del:mi, where: dist=distance bewteeen each complex varant; len=length of each complex variant; s=number of snps; ins=number of insertions; del=number of deletions; mi=max indel length",
        metavar="LIST1[,LIST2,...]"
    )

    subparser_mutate_fasta.add_argument(
        "fasta_in",
        help="FASTA filename of genome to be mutated",
    )

    subparser_mutate_fasta.add_argument(
        "outprefix",
        help="Prefix of output files",
    )

    subparser_mutate_fasta.set_defaults(func=simutator.tasks.mutate_fasta.run)

    args = parser.parse_args()
    logging.basicConfig(
        format=f"[%(asctime)s simutator %(levelname)s] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    log = logging.getLogger()
    if args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
