#!/usr/bin/env python
import sys
from pathlib import Path
import shutil
import argparse
import traceback
import re
import yaml
import pandas as pd, numpy as np
from dataclasses import dataclass, field
import csv

mapping_types = {
    "pairID": str,
    "fwdprimer_name": str,
    "fwdprimer": str,
    "revprimer_name": str,
    "revprimer": str,
    "prior_size": "Int16",
    "amp_length_primers": str,
    "amplicon_length_total": str,
    "max_length_total": "Int16",
    "max_length_primers": "Int16",
    "variable_length": str,
    "trimFor": "Int16",
    "trimRev": "Int16",
    "truncFor": "Int16",
    "truncRev": "Int16",
    "maxLen": "Int16",
    "minLen": "Int16",
    "maxEEFor": "Int16",
    "maxEERev": "Int16",
    "minOverlap": "Int16",
    "trimOverhang": str,
    "justConcatenate": str,
    "minMergedLength": "Int16",
    "maxMergedLength": "Int16",
    "reference": str,
    "species": str,
    "removeBimeraDenovoOptions": str,
    "runTree": str,
    "taxassignment": str,
    "notes": str,
}


# TODO: decide on best decorator for this, maybe @dataclass?
@dataclass
class AmpliconSet(object):
    """
    Simple data class to hold amplicon data along with
    associated functions for guessing settings given a
    specific read length.
    """

    pairID: str
    forward_primer_name: str
    reverse_primer_name: str
    forward_primer: str
    reverse_primer: str
    max_length_primers: int = None
    variable_length: bool = False
    notes: str = None


# TODO: set up a minimal parameters object for DADA2 runs
# TODO: set global defaults for DADA2 in this
@dataclass
class TADAParams(object):
    """
    TADA Params (in flux, beware!)
    """

    _amplicon: AmpliconSet = None
    _read_length: int = 250
    _min_allowed_overlap: int = 10
    _paired_end: bool = True  # at the moment this is all we support
    _strategy: str = "paired"
    project: str = None
    trimFor: int = None
    trimRev: int = None
    fwdprimer: str = None
    revprimer: str = None
    pool: str = "pseudo"
    email: str = None
    truncFor: int = None
    truncRev: int = None
    input: str = None
    skip_dadaQC: bool = False
    skip_multiQC: bool = False
    precheck: bool = True
    check_merging: bool = False
    outdir: str = None
    maxEEFor: int = 2  # DADA2 default = Inf
    maxEERev: int = 2  # DADA2 default = Inf
    truncQ: int = 2  # DADA2 default = 2
    maxLen: int = None
    minLen: int = None
    minOverlap: int = None
    maxMismatch: int = None
    trimOverhang: bool = None
    justConcatenate: bool = False
    removeBimeraDenovoOptions: str = None
    qualityBinning: bool = False
    minMergedLength: int = 50
    maxMergedLength: int = None
    reference: str = None
    taxassignment: str = None
    aligner: str = "DECIPHER"
    runTree: str = None
    species: str = None
    toQIIME2: bool = True
    idType: str = "md5"
    interactiveMultiQC: bool = True

    def __post_init__(self):
        """
        This is kinda clunky, but the general idea is how
        we handle certain settings not in the parameters
        """
        if pd.isna(self.trimFor) and pd.isna(self.trimRev):
            self.set_trimming()
        if pd.isna(self.truncFor) and pd.isna(self.truncRev):
            self.guess_trunc()
        if pd.isna(self.minOverlap):
            self.guess_minOverlap()
        self.check_overhang()
        if pd.isna(self.outdir):
            self.outdir = self._amplicon.pairID
        if pd.isna(self.project):
            self.project = self._amplicon.pairID

    def set_trimming(self):
        """
        If this isn't set, default to the instance
        """
        self.fwdprimer = self._amplicon.forward_primer
        self.revprimer = self._amplicon.reverse_primer
        self.trimFor = len(self._amplicon.forward_primer)
        self.trimRev = len(self._amplicon.reverse_primer)

    def guess_trunc(self):
        """
        Guess truncation parameters based on AmpliconSet
        """
        # =IF(M3<'PrimerDB Settings'!$B$1, M3, 0)
        truncate = (
            self._amplicon.max_length_primers
            if (
                self._amplicon.max_length_primers
                and self._amplicon.max_length_primers < self._read_length
            )
            else 0
        )
        self.truncFor = truncate
        if self._paired_end:
            self.truncRev = truncate

    def guess_minOverlap(self):
        """
        Guess overlap based on AmpliconSet
        """
        forlen = self.truncFor if self.truncFor > 0 else self._read_length
        revlen = self.truncRev if self.truncRev > 0 else self._read_length
        if self._amplicon.max_length_primers and self._amplicon.max_length_primers <= (
            forlen + revlen
        ):
            # the fudge factor here is completely arbitrary, but it's to allow
            # for some level of variance in the overlapped region
            minovl = int(
                (forlen + revlen)  # total of truncated or full-length reads
                - self._amplicon.max_length_primers  # amplicon w/ primers
                - (self._read_length * 0.1)  # fudge factor
            )
            if self._amplicon.variable_length == "True" and minovl < 50:
                self._strategy = "variable"
            # we don't want to go less than 5-10 if possible, and no more than
            # 150 (this may be something we parameterize)
            if minovl <= 10:
                minovl = 10

            if minovl >= 150:
                minovl = 150

            self.minOverlap = minovl
        else:
            # no overlap detected but maybe the amplicon length isn't known
            self._strategy = "single" if self._amplicon.max_length_primers else "check"

    def check_overhang(self):
        """
        Check if there is a potential overhang. This pops up if we manually
        set truncFor/Rev (it's defined in the sample sheet)
        """
        self.trimOverhang = (
            self._amplicon.max_length_primers
            and self._amplicon.max_length_primers <= (self.truncFor + self.truncRev)
        )


def main():
    parser = argparse.ArgumentParser(
        description="""
    Simple script to take data from a primer-sorted Fluidigm run and generate:

    1) formatted samplesheet with sample IDs and the full path to the FASTQ
       data
    2) a general configuration file for all the runs (can be passed in from
       the command line), and (optionally but recommended)

    3) a parameters file using information from a spreadsheet of primer
    combinations.

    This assumes the directory structure the '-s' option points to has a folder
    with sequences split into subfolders by primer pair; this primer pair is
    used for optionally identifying from the mapping file which parameters to
    use.

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
        "--single_end",
        action="store_true",
        help="""
        Common config file for runs (leave blank to use baseline config)
        """,
    )

    parser.add_argument(
        "--quality_binning",
        action="store_true",
        help="""
        Common config file for runs (leave blank to use baseline config)
        """,
    )

    args = parser.parse_args()

    print("Step 1: get sequence files")
    seq_files = parse_seqdata(args.seqdata)

    print("Step 2: get primer pair data")
    mapping_file = None

    if args.fluidigm:
        mapping_file = parse_fluidigm_mapping(
            mapping=args.fluidigm, email=args.email, binning=args.quality_binning
        )

    print("Step 3: set up workspace")
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


def parse_fluidigm_mapping(mapping=None, email=None, binning=False):
    """
    Parsing the mapping files of attributes to primers;
    very little validation at the moment
    """
    fl_mapping = {}

    # the current version of the 'database' is an Excel table
    fluidigm_data = pd.read_excel(
        mapping,
        index_col=0,
        dtype=mapping_types,
        comment="#",
        sheet_name="PrimerDB",
    )

    fluidigm_data = fluidigm_data.replace({np.nan: None})
    # TODO: switch to vectorization (using this to debug)

    for idx, row in fluidigm_data.iterrows():
        tada_params = row_to_TADAParams(row, email=email)
        if binning:
            tada_params.qualityBinning = True
        fl_mapping[tada_params._amplicon.pairID] = tada_params
    return fl_mapping


def mk_double_quote(dumper, data):
    dumper.represent_scalar("tag:yaml.org,2002:str", data, style='"')


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
        if pair not in mapping_file:
            print("Pair %s doesn't match anything in the database?" % pair)
            continue
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

        # parameters
        parampath = pairpath.joinpath("params.%s.yml" % pair)
        with open(parampath, "w") as param_yaml:
            mapping = {}
            if mapping_file and pair in mapping_file:
                mapping_file[pair].input = str(samplesheet)
                # testing out filtering private and null/NA
                mapping = mapping_file[pair].__dict__
                for k, v in list(mapping.items()):
                    if k.startswith("_") or v is None:
                        del mapping[k]

                # TODO: we can set up strategy-specific settings here
                # if mapping_file[pair]._strategy != "paired":
            else:
                print(f"Pair {pair} not found, using stub")
                mapping = {p: None for p in params}

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
