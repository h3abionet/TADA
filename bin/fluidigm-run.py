#!/usr/bin/env python
import sys
from pathlib import Path
import shutil
import argparse
import traceback
import re
import yaml
import pandas as pd
from dataclasses import dataclass

# import csv

# import itertools


# TODO: decide on best decorator for this
@dataclass
class AmpliconSet(object):
    """
    Simple data class to hold amplicon data along with
    associated functions for guessing settings given a
    specific read length.

    General primer info:
        paired_end=True,
        forward_primer=None,
        reverse_primer=None,
        amplicon_length=None,
        read_length=[0,0]
    Trimming behavior:
        trunc=[0,0],
        maxEE=[2,2],
    Filtering (trimming):
        maxLength=0,
        minLength=0,
        minOverlap=0,
        trimOverhang=False,
        justConcatenate=False,
        minMergedLength=50,
        maxMergedLength=None

    """

    paired_end: bool = True
    forward_primer: str = None
    reverse_primer: str = None
    amplicon_length: int = None


# TODO: set up a minimal parameters object for DADA2 runs
# TODO: set global defaults for DADA2 in this
class DadaParams(object):
    """docstring for ClassName"""

    def __init__(self, arg):
        super(ClassName, self).__init__()
        self.arg = arg


def main():
    parser = argparse.ArgumentParser(
        description="""
    Simple script to take data from a primer-sorted Fluidigm run and generate:
    1) formatted samplesheet with sample IDs and the full path to the FASTQ data
    2) a general configuration file for all the runs (can be passed in from 
       the command line)

    and (optionally but recommended)

    3) a parameters file using information from a spreadsheet of 
    primer combinations

    This assumes the directory structure the '-s' option points to has a folder
    with sequences split into subfolders by primer pair; this primer pair is used
    for optionally identifying from the mapping file which parameters to use.

    Note that some will be predictive values if defaults are not present.
    """
    )
    parser.add_argument(
        "-f",
        "--fluidigm",
        required=False,
        help="""
        Fluidigm mapping file
        """,
        default=None,
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
        "-e",
        "--email",
        default=None,
        help="""
        Email address for reports
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
        "-g",
        "--guess",
        required=False,
        help="""
        If using parameter mapping and some fields are blank, 
        guess settings (use with caution)
        """,
        default=None,
    )
    parser.add_argument(
        "--paired",
        default=True,
        help="""
        Paired reads (set to False if single-end
        """,
    )
    args = parser.parse_args()

    # print("Step 1: get sequence files")
    seq_files = parse_seqdata(args.seqdata)

    # print("Step 2: get primer pair data")
    mapping_file = None
    # if args.fluidigm:
    #     mapping_file = parse_fluidigm_mapping(args.fluidigm)

    # print("Step 3: set up workspace")
    setup_workspace(args, seq_files, mapping_file)


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
        read_match = re.match(f"{pair_name}-(\S+)_[ATGC]+_(R\d)", sample_file)
        if read_match:
            sn = read_match.group(1)
            if sn not in seq_files[pair_name]:
                seq_files[pair_name][sn] = {}
            seq_files[pair_name][sn][read_match.group(2)] = str(f)
        else:
            print("Unmatched file: ", sample_file)
            # try to get the name of the sample
    return seq_files


def parse_fluidigm_mapping(mapping):
    """
    Parsing the mapping files of attributes to primers;
    very little validation at the moment
    """
    fl_mapping = {}
    # to_bool = ["justConcatenate", "trimOverhang", "variable"]
    # to_int = [
    #     "maxEEFor",
    #     "maxEERev",
    #     "minOverlap",
    #     "truncFor",
    #     "truncRev",
    #     "trimFor",
    #     "trimRev",
    # ]

    # the current version of the 'database' is an Excel table
    fluidigm_data = pd.read_excel(mapping, index_col=1, sheet_name="PrimerDB")

    print(fluidigm_data[0:3])
    # csv_as_string = xlsx_data.to_csv(index=False)
    # reader = csv.DictReader(csv_as_string.splitlines())
    # fl_reader = csv.DictReader(
    #     filter(lambda row: row[0] != "#", csv_as_string.splitlines()),
    # )
    # for row in fl_reader:
    #     # we want simple lookup table 'PrimerPairID'
    #     tmp = dict(itertools.islice(row.items(), 2, None))
    #     for i in to_bool:
    #         tmp[i] = bool(tmp[i])
    #     for i in to_int:
    #         if tmp[i] == "":
    #             tmp[i] = 0
    #         tmp[i] = int(tmp[i])
    #     fl_mapping[row["PrimerPairID"]] = tmp
    # return fl_mapping


global dict_keys


def mk_double_quote(dumper, data):
    dumper.represent_scalar("tag:yaml.org,2002:str", data, style='"')


# def setup_workspace(args, seq_files, mapping_file):


def setup_workspace(args, seq_files, mapping_file):
    """
    Iterate through the seq files by primer pair
    and check them against mapping file
    """
    workdir = Path(args.workdir).resolve()
    params = [
        "trimFor",
        "trimRev",
        "fwdprimer",
        "revprimer",
        "truncFor",
        "truncRev",
        "minOverlap",
        "maxEEFor",
        "maxEERev",
        "justConcatenate",
        "trimOverhang",
        "reference",
        "species",
        "variable",
    ]
    for pair in seq_files.keys():
        # if pair not in mapping_file:
        #     print("Pair %s doesn't match anything in the database?" % pair)
        #     continue
        pairpath = workdir.joinpath(pair)

        # Generate primer pair workspace
        try:
            pairpath.mkdir(parents=True)
        except FileExistsError as e:
            raise e
        except FileNotFoundError as e:
            raise e

        # set up sample sheet
        # TODO: support explicit or implicit single-end
        samplesheet = pairpath.joinpath("samplesheet.%s.csv" % pair)
        with open(samplesheet, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["sample", "fastq_1", "fastq_2"])
            for sample in seq_files[pair].keys():
                writer.writerow(
                    [
                        sample,
                        seq_files[pair][sample]["R1"],
                        seq_files[pair][sample]["R2"],
                    ]
                )

        # copy config file if available
        configpath = Path(args.config)
        if configpath.is_file():
            print("Found file, copying to ", str(pairpath))
            shutil.copyfile(configpath, pairpath.joinpath(configpath.name))
        else:
            print("Skipping config file")

        # parameters
        parampath = pairpath.joinpath("params.%s.yml" % pair)
        with open(parampath, "w") as param_yaml:
            # if pair == "V4_515F_New_V4_806R_New":
            #     # add project name, email, input (samplesheet), name,
            #     # outdir, skip_dadaQC, interactiveMultiQC=FALSE, pool,
            #     # aligner, runTree, toQIIME2, idType, qualityBinning
            opts = {
                "pool": True,
                "project": pair,
                "name": pair,
                "outdir": pair,
                "skip_dadaQC": True,
                "skip_multiQC": True,
                "interactiveMultiQC": False,
                "pool": "pseudo",
                "aligner": "DECIPHER",
                "runTree": "fasttree",
                "toQIIME2": True,
                "idType": "md5",
                "qualityBinning": True,
                "email": args.email,
                "input": str(samplesheet),
            }

            mapping = {}
            if mapping_file:
                if pair in mapping_file:
                    for p in params:
                        mapping = {p: mapping_file[pair][p] for p in params}
            else:
                mapping = {p: None for p in params}

            mapping.update(opts)

            yaml.dump(mapping, param_yaml)


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
