#!/usr/bin/env python
import sys
from pathlib import Path
import shutil
import argparse
import traceback
import re
from dataclasses import dataclass, field
import csv


def main():
    parser = argparse.ArgumentParser(
        description="""
    Simple script to take data from a primer-sorted Fluidigm run and generate
    formatted samplesheets with sample IDs and the full path to the FASTQ data
    
    This assumes the directory structure the '-s' option points to has a folder
    with sequences split into subfolders named by primer pair.
    """
    )

    parser.add_argument(
        "-s",
        "--seqdata",
        required=True,
        help="""
        Location of sequence data
        """,
    )
    parser.add_argument(
        "-w",
        "--workdir",
        required=True,
        help="""
        Location of directory to set up jobs
        """,
    )
    parser.add_argument(
        "-c",
        "--config",
        default="",
        help="""
        Common config file for runs (leave blank to use baseline config)
        """,
    )
    parser.add_argument(
        "--single_end",
        action="store_true",
        help="""
        If single-end samplesheet is needed
        """,
    )
    args = parser.parse_args()

    print("Step 1: get sequence files")
    seq_files = parse_seqdata(args.seqdata)

    print("Step 2: set up workspace")
    setup_workspace(args, seq_files)


def parse_seqdata(dir):
    """
    Parse a given directory for folders (named by primer pair)
    and files for each read.  Note this could also be used for
    R1 only
    """
    # structure: primer pair - sample name - read pair
    seq_files = {}
    # we want absolute paths and only files in the top directories
    seqdir = Path(dir).resolve().glob("*/*.fastq.gz")
    for f in seqdir:
        pair_name = f.parts[-2]
        if pair_name not in seq_files:
            seq_files[pair_name] = {}
        sample_file = f.parts[-1]
        read_match = re.match(f"{pair_name}-(\\S+)_[ATGC]+_(R\\d)", sample_file)
        if read_match:
            sn = read_match.group(1)
            if sn not in seq_files[pair_name]:
                seq_files[pair_name][sn] = {}
            seq_files[pair_name][sn][read_match.group(2)] = str(f)
        else:
            print("Unmatched file: ", sample_file)
            # try to get the name of the sample
    return seq_files


def row_to_TADAParams(row, email=None):
    # scrub inputs to remove anything NA or None
    amp_data = row[
        [
            "pairID",
            "forward_primer",
            "forward_primer_name",
            "reverse_primer",
            "reverse_primer_name",
            "max_length_primers",
            "variable_length",
            "notes",
        ]
    ]

    # empty cells clobber the instance defaults, so we scrub them
    for k, v in list(amp_data.items()):
        if v is None:
            del amp_data[k]
    ampset = AmpliconSet(**amp_data)

    params_data = row[
        [
            "truncFor",
            "truncRev",
            "maxEEFor",
            "maxEERev",
            "maxLen",
            "minLen",
            "minOverlap",
            "trimOverhang",
            "justConcatenate",
            "minMergedLength",
            "maxMergedLength",
            "reference",
            "runTree",
            "taxassignment",
            "species",
        ]
    ]
    # empty cells clobber the instance defaults, so we scrub them
    for k, v in list(params_data.items()):
        if v is None:
            del params_data[k]
    tada_params = TADAParams(_amplicon=ampset, email=email, **params_data)
    return tada_params


def mk_double_quote(dumper, data):
    dumper.represent_scalar("tag:yaml.org,2002:str", data, style='"')


def setup_workspace(args, seq_files):
    """
    Iterate through the seq files by primer pair
    and check them against mapping file
    """
    workdir = Path(args.workdir).resolve()

    for pair in seq_files.keys():
        # if pair not in mapping_file:
        #     print("Pair %s doesn't match anything in the database?" % pair)
        #     continue
        pairpath = workdir.joinpath(pair)

        # Generate primer pair workspace
        try:
            # TODO: default to not overwriting, add bool flag
            pairpath.mkdir(parents=True, exist_ok=True)
        except FileExistsError as e:
            raise e
        except FileNotFoundError as e:
            raise e

        # set up sample sheet
        # TODO: support explicit or implicit single-end
        samplesheet = pairpath.joinpath("samplesheet.%s.pe.csv" % pair)
        with open(samplesheet, "w", newline="") as csvfile:
            writer = csv.writer(csvfile, lineterminator="\n")
            writer.writerow(["sample", "fastq_1", "fastq_2"])
            for sample in seq_files[pair].keys():
                writer.writerow(
                    [
                        sample,
                        seq_files[pair][sample]["R1"],
                        seq_files[pair][sample]["R2"],
                    ]
                )

        if args.single_end:
            print("Generating single-end sample sheet just in case...")
            samplesheet_se = pairpath.joinpath("samplesheet.%s.se.csv" % pair)
            with open(samplesheet_se, "w", newline="") as se_csvfile:
                writer = csv.writer(se_csvfile, lineterminator="\n")
                writer.writerow(["sample", "fastq_1", "fastq_2"])
                for sample in seq_files[pair].keys():
                    writer.writerow(
                        [
                            sample,
                            seq_files[pair][sample]["R1"],
                            "",
                        ]
                    )

        # copy config file if available
        configpath = Path(args.config)
        if configpath.is_file():
            print("Found file, copying to ", str(pairpath))
            shutil.copyfile(configpath, pairpath.joinpath(configpath.name))
        else:
            print("Skipping config file")


if __name__ == "__main__":
    try:
        main()
        sys.exit(0)
    except KeyboardInterrupt as e:
        raise e
    except SystemExit as e:
        raise e
    except Exception as e:
        print(e)
        traceback.print_exc()
        sys.exit(1)
