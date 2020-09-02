#!/usr/bin/env python
# coding: utf-8

import argparse
import sys
from datetime import datetime

from logzero import logger

#from . import util
from . import misc
from . import cnv_arithmetics as ca
from . import matchscores
from . import cnv_calling
from . import visualize_scores
from . import workflow
from . import __version__


def get_parser():
    """Return argparse command line parser."""
    parser = argparse.ArgumentParser(
        description="clearCNV can compute matchscores, CNV-calls and visualizations.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers()

    # =========================================================================
    #  matchscores
    # =========================================================================
    parser_matchscores = subparsers.add_parser(
        "matchscores", description="Matchscore calculation script."
    )
    parser_matchscores.add_argument(
        "-p", "--panel", help="Name of the data set (or panel)", required=True, type=str
    )
    parser_matchscores.add_argument(
        "-c", "--coverages", help="Coverages file in tsv format", required=True, type=str
    )
    parser_matchscores.add_argument(
        "-m", "--matchscores", help="Output matchscores.tsv file", required=True, type=str
    )
    parser_matchscores.add_argument(
        "-x",
        "--expected_artefacts",
        help="Expected ratio of CNVs or artefacs in target fragment counts",
        required=False,
        type=float,
        default=0.02,
    )
    parser_matchscores.add_argument(
        "--cores",
        help="Number of cpu cores used in parallel processing. Default: determined automatically.",
        required=False,
        type=int,
        default=0,
    )
    parser_matchscores.set_defaults(func=matchscores.matchscores)

    # =========================================================================
    #  cnv calling
    # =========================================================================
    parser_cnv_calling = subparsers.add_parser(
        "cnv_calling",
        description=(
            "CNV calling script. Output is a single file in tsv format containing a list of CNV "
            "calls sorted by score. Some quality control plots are added to the analysis "
            "directory in the process."
        ),
    )
    parser_cnv_calling.add_argument(
        "-p", "--panel", help="Name of the data set(or panel)", required=True, type=str
    )
    parser_cnv_calling.add_argument(
        "-c", "--coverages", help="Coverages file in tsv format", required=True, type=str
    )
    parser_cnv_calling.add_argument(
        "-a",
        "--analysis_directory",
        help="Path to the directory, where analysis files are stored",
        required=True,
        type=str,
    )
    parser_cnv_calling.add_argument(
        "-m",
        "--matchscores",
        help="matchscores.tsv file generated with matchscores.py",
        required=True,
        type=str,
    )
    parser_cnv_calling.add_argument(
        "-C",
        "--cnv_calls",
        help="Output cnv.calls file formatted in tsv format",
        required=True,
        type=str,
    )
    parser_cnv_calling.add_argument(
        "-r",
        "--ratio_scores",
        help="Output ratio scores file in tsv format. Best kept together with cnv.calls",
        required=True,
        type=str,
    )
    parser_cnv_calling.add_argument(
        "-z",
        "--z_scores",
        help="Output z-scores file in tsv format. Best kept together with cnv.calls",
        required=True,
        type=str,
    )
    parser_cnv_calling.add_argument(
        "-x",
        "--expected_artefacts",
        help="Expected ratio of CNVs or artefacs in target fragment counts",
        required=False,
        type=float,
        default=0.02,
    )
    parser_cnv_calling.add_argument(
        "-u",
        "--minimum_sample_score",
        help="A lower threshold results in better fitting, but smaller calling groups",
        required=False,
        type=float,
        default=0.15,
    )
    parser_cnv_calling.add_argument(
        "-g",
        "--minimum_group_sizes",
        help="Minimum group size per CNV calling group per match scores",
        required=False,
        type=int,
        default=33,
    )
    parser_cnv_calling.add_argument(
        "-s",
        "--sensitivity",
        help="A higher sensitivity results in more CNV calls. Can only be 0.0 <= sens <= 1.0",
        required=False,
        type=float,
        default=0.7,
    )
    parser_cnv_calling.add_argument(
        "--cores",
        help="Number of cpu cores used in parallel processing. Default: determined automatically.",
        required=False,
        type=int,
        default=0,
    )
    parser_cnv_calling.set_defaults(func=cnv_calling.cnv_calling)

    # =========================================================================
    #  visualize
    # =========================================================================
    parser_visualize = subparsers.add_parser(
        "visualize",
        description=(
            "The visualization script creates html files containing heatmap-like matrices "
            "containing the ratios aligned with mappability, GC-content and target size so "
            "that CNVs can be visually identified and evaluated easily."
        ),
    )
    parser_visualize.add_argument(
        "-a",
        "--analysis_directory",
        help="Path to the directory, where analysis files are stored",
        required=True,
        type=str,
    )
    parser_visualize.add_argument(
        "-r",
        "--ratio_scores",
        help="Ratio scores file in tsv format, generated in cnv_calling.py",
        required=True,
        type=str,
    )
    parser_visualize.add_argument(
        "-z",
        "--z_scores",
        help="Z-scores file in tsv format, generated in cnv_calling.py",
        required=True,
        type=str,
    )
    parser_visualize.add_argument(
        "-n",
        "--annotated",
        help="CALL_GROUPS.py analysis_directory z_scores.tsv ratio_scores.tsv annotations.bed",
        required=True,
        type=str,
    )
    parser_visualize.add_argument(
        "-s",
        "--size",
        help="Rough number of targets in each visualization. Defaults at 1000.",
        required=False,
        type=int,
        default=1000,
    )
    parser_visualize.set_defaults(func=visualize_scores.visualize)

    # =========================================================================
    #  merge bed
    # =========================================================================
    parser_merge_bed = subparsers.add_parser(
        "merge_bed",
        description=("Merges bed files to non-overlapping intervals."),
    )
    parser_merge_bed.add_argument(
        "-i",
        "--infile",
        help="Path to the original .bed file.",
        required=True,
        type=str,
    )
    parser_merge_bed.add_argument(
        "-o",
        "--outfile",
        help="Output path to the merged .bed file.",
        required=True,
        type=str,
    )
    parser_merge_bed.set_defaults(func=misc.merge_bedfile)

    # =========================================================================
    #  coverage
    # =========================================================================
    parser_rtbeds = subparsers.add_parser(
        "coverage",
        description=(
            "wrapper for bedtools multicov. Creates an .rtbed file which contains the read depth coverage per target."
        ),
    )
    parser_rtbeds.add_argument(
        "-i",
        "--inbam",
        help="Path to the .bam file of the sample.",
        required=True,
        type=str,
    )
    parser_rtbeds.add_argument(
        "-b",
        "--bedfile",
        help="Path to the merged .bed file.",
        required=True,
        type=str,
    )
    parser_rtbeds.add_argument(
        "-r",
        "--rtbed",
        help="Output file in .rtbed format.",
        required=True,
        type=str,
    )
    parser_rtbeds.set_defaults(func=misc.count_pair)

    # -------------------------------------------------------------------------
    # merge rtbeds

    parser_merge_rtbeds = subparsers.add_parser(
        "merge_coverages",
        description=("Merges all .rtbed files into one table in .tsv format."),
    )
    parser_merge_rtbeds.add_argument(
        "-b",
        "--bedfile",
        help="Path to the merged .bed file.",
        required=True,
        type=str,
    )
    parser_merge_rtbeds.add_argument(
        "-r",
        "--rtbeds",
        help="Input .rtbed file paths.",
        required=True,
        nargs="+",
        type=str,
    )
    parser_merge_rtbeds.add_argument(
        "-c",
        "--coverages",
        help="Output table in .tsv format.",
        required=True,
        type=str,
    )
    parser_merge_rtbeds.set_defaults(func=misc.merge_rtbeds)

    # =========================================================================
    #  annotations
    # =========================================================================
    parser_annotations = subparsers.add_parser(
        "annotations",
        description=("Creates annotations file."),
    )
    parser_annotations.add_argument(
        "-r",
        "--reference",
        help="Path to the genomic reference.",
        required=True,
        type=str,
    )
    parser_annotations.add_argument(
        "-b",
        "--bedfile",
        help="Path to the merged .bed file.",
        required=True,
        type=str,
    )
    parser_annotations.add_argument(
        "-k",
        "--kmer_align",
        help="Path to aligned k-mers (mappability) file in .bed format.",
        required=True,
        type=str,
    )
    parser_annotations.add_argument(
        "-a",
        "--annotations",
        help="Output file in .bed format.",
        required=True,
        type=str,
    )
    parser_annotations.set_defaults(func=misc.create_annotations_file)

    # =========================================================================
    #  workflow
    # =========================================================================
    parser_workflow = subparsers.add_parser(
        "workflow",
        description=("Complete workflow that runs snakemake internally."),
    )
    parser_workflow.add_argument(
        "-w",
        "--workdir",
        help="Path to the snakemake workdir.",
        required=True,
        type=str,
    )
    parser_workflow.add_argument(
        "-p",
        "--panelname",
        help="name of the panel or dataset.",
        required=True,
        type=str,
    )
    parser_workflow.add_argument(
        "-r",
        "--reference",
        help="Path to the genomic reference.",
        required=True,
        type=str,
    )
    parser_workflow.add_argument(
        "-b",
        "--bamsfile",
        help="Path to a .txt file that contains all paths to all used .bam files.",
        required=True,
        type=str,
    )
    parser_workflow.add_argument(
        "-d",
        "--bedfile",
        help="Path to the .bed file.",
        required=True,
        type=str,
    )
    parser_workflow.add_argument(
        "-k",
        "--kmer_align",
        help="Path to aligned k-mers (mappability) file in .bed format.",
        required=True,
        type=str,
    )
    parser_workflow.set_defaults(func=workflow.workflow)


    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
    else:
        logger.info("Starting analysis...")
        args.func(args)
        logger.info("All done. Have a nice day!")


if __name__ == "__main__":
    main()
