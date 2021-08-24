#!/usr/bin/env python

import argparse
import os
import tempfile

from logzero import logger

from . import misc
from . import matchscores
from . import cnv_calling
from . import visualize_scores
from . import workflow_cnv_calling
from . import workflow_reassignment
from . import __version__

#: The executables required for running clear-CNV.
REQUIRED_EXECUTABLES = ("bedops", "bedtools", "sort")


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
    parser_matchscores.add_argument(
        "--fast",
        help="If set, clearCNV will speed up matchscore calculation by taking only at most 2.000 targets per sample into account.",
        required=False,
        dest="fast",
        action="store_true",
    )
    parser_matchscores.set_defaults(fast=False, func=matchscores.matchscores)

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
        help="Output cnv_calls.tsv file formatted in tsv format",
        required=True,
        type=str,
    )
    parser_cnv_calling.add_argument(
        "-r",
        "--ratio_scores",
        help="Output ratio scores file in tsv format. Best kept together with cnv_calls.tsv",
        required=True,
        type=str,
    )
    parser_cnv_calling.add_argument(
        "-z",
        "--z_scores",
        help="Output z-scores file in tsv format. Best kept together with cnv_calls.tsv",
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
        "--sample_score_factor",
        help="The factor u multiplied with the median sample score to define sample groups. u should range between 1.0 < u < 5.0. Default is 2.0.",
        required=False,
        type=float,
        default=2.0,
    )
    parser_cnv_calling.add_argument(
        "-g",
        "--minimum_group_sizes",
        help="Group size per CNV calling group per match scores. Default is 30.",
        required=False,
        type=int,
        default=30,
    )
    parser_cnv_calling.add_argument(
        "-s",
        "--zscale",
        help="A higher z-scale results in more CNV calls. Should only be 0.0 <= zscale <= 2.0. Default is 0.65.",
        required=False,
        type=float,
        default=0.65,
    )
    parser_cnv_calling.add_argument(
        "--del_cutoff",
        help="A hard threshold on the ratio score for deletions. Default is 0.75.",
        required=False,
        type=float,
        default=0.75,
    )
    parser_cnv_calling.add_argument(
        "--dup_cutoff",
        help="A hard threshold on the ratio score for duplications. Default is 1.35.",
        required=False,
        type=float,
        default=1.35,
    )
    parser_cnv_calling.add_argument(
        "--trans_prob",
        help="Transition probability of the HMM to change state from WT to CNV. Default is 0.001. Lower values prefer longer CNVs, higher values prefer shorter CNVs.",
        required=False,
        type=float,
        default=0.001,
    )
    parser_cnv_calling.add_argument(
        "--plot_regions",
        dest="plot_regions",
        help="If set, the CNV calling script plots heatmaps of the CNV called regions with the corresponding sample group taken as background.",
        required=False,
        action="store_true",
    )
    parser_cnv_calling.set_defaults(plot_regions=False)
    # parser_cnv_calling.add_argument(
    #    "--plot_regions",
    #    help="If true, the CNV calling script plots heatmaps of the regions.",
    #    required=False,
    #    type=bool,
    #    default=False,
    # )
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
        help="BED-formatted annotatins file generated with clearCNV annotations. E.g. annotations.bed.",
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
    #  reassignment
    # =========================================================================

    # prepare - maybe exclude in future
    parser_prepare_reassign = subparsers.add_parser(
        "prepare_reassignment",
        description=(
            "Prepares the necessary files to perform panel reassignment on a set of .bed files and corresponding .bam file data sets."
        ),
    )
    parser_prepare_reassign.add_argument(
        "-m",
        "--metafile",
        help="Path to the file containing the meta information. It is a .tsv of the scheme -panel bams.txt bed.bed. \
            It aligns each desired panel name with the corresponding .bam files and the .bed file.",
        required=True,
        type=str,
    )
    parser_prepare_reassign.add_argument(
        "-b",
        "--bamsfile",
        help="Output .txt file. It will contain the distinct set of all given .bam file paths.",
        required=True,
        type=str,
    )
    parser_prepare_reassign.add_argument(
        "-d",
        "--bedfile",
        help="Output .bed file. It will contain the merged union of all given .bed files.",
        required=True,
        type=str,
    )
    parser_prepare_reassign.set_defaults(func=misc.prepare_reassignment)

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
    #  workflow CNV calling
    # =========================================================================
    parser_workflow_cnv_calling = subparsers.add_parser(
        "workflow_cnv_calling",
        description=(
            "Complete CNV-calling snakemake-workflow to run on data from a single sequencing panel (bed-file)."
        ),
    )
    parser_workflow_cnv_calling.add_argument(
        "-w",
        "--workdir",
        help="Path to the snakemake workdir.",
        required=True,
        type=str,
    )
    parser_workflow_cnv_calling.add_argument(
        "-p",
        "--panelname",
        help="name of the panel or dataset.",
        required=True,
        type=str,
    )
    parser_workflow_cnv_calling.add_argument(
        "-r",
        "--reference",
        help="Path to the genomic reference.",
        required=True,
        type=str,
    )
    parser_workflow_cnv_calling.add_argument(
        "-b",
        "--bamsfile",
        help="Path to a .txt file that contains all paths to all used .bam files.",
        required=True,
        type=str,
    )
    parser_workflow_cnv_calling.add_argument(
        "-d",
        "--bedfile",
        help="Path to the BED file that contains all targets of the sequencing panel. It must contain the columns: chr, start, end, gene.",
        required=True,
        type=str,
    )
    parser_workflow_cnv_calling.add_argument(
        "-k",
        "--kmer_align",
        help="Path to aligned k-mers (mappability) file in .bed format.",
        required=True,
        type=str,
    )
    parser_workflow_cnv_calling.add_argument(
        "-c",
        "--cores",
        help="Number of used cores in snakemake workflow. Default is 32.",
        required=False,
        type=int,
        default=32,
    )
    # options
    parser_workflow_cnv_calling.add_argument(
        "--expected_artefacts",
        help="Expected ratio of CNVs or artefacs in target fragment counts. Default is 0.02",
        required=False,
        type=float,
        default=0.02,
    )
    parser_workflow_cnv_calling.add_argument(
        "--sample_score_factor",
        help="The factor u multiplied with the median sample score to define sample groups. u should range between 1.0 < u < 5.0. Default is 2.0.",
        required=False,
        type=float,
        default=2.0,
    )
    parser_workflow_cnv_calling.add_argument(
        "--minimum_group_sizes",
        help="Minimum group size per CNV calling group per match scores. Default is 20.",
        required=False,
        type=int,
        default=20,
    )
    parser_workflow_cnv_calling.add_argument(
        "--zscale",
        help="A higher z-scale results in more CNV calls. Should only be 0.0 <= zscale <= 2.0. Default is 0.65.",
        required=False,
        type=float,
        default=0.65,
    )
    parser_workflow_cnv_calling.add_argument(
        "--size",
        help="Rough number of targets in each visualization. Default is 2000.",
        required=False,
        type=int,
        default=2000,
    )
    parser_workflow_cnv_calling.add_argument(
        "--del_cutoff",
        help="A hard threshold on the ratio score for deletions. Default is 0.75.",
        required=False,
        type=float,
        default=0.75,
    )
    parser_workflow_cnv_calling.add_argument(
        "--dup_cutoff",
        help="A hard threshold on the ratio score for duplications. Default is 1.35.",
        required=False,
        type=float,
        default=1.35,
    )
    parser_workflow_cnv_calling.add_argument(
        "--trans_prob",
        help="Transition probability of the HMM to change state from WT to CNV. Default is 0.001. Lower values prefer longer CNVs, higher values prefer shorter CNVs.",
        required=False,
        type=float,
        default=0.001,
    )
    parser_workflow_cnv_calling.add_argument(
        "--plot_regions",
        dest="plot_regions",
        help="If set, the CNV calling script plots heatmaps of the CNV called regions with the corresponding sample group taken as background.",
        required=False,
        action="store_true",
    )
    parser_workflow_cnv_calling.set_defaults(plot_regions=False)
    parser_workflow_cnv_calling.add_argument(
        "--cluster_configfile",
        help="Path to the cluster config file.",
        required=False,
        type=str,
        default="",
    )

    parser_workflow_cnv_calling.add_argument(
        "--drmaa_mem",
        help="Number of megabytes used with drmaa. Suggested is 16000. If --drmaa_mem and --drmaa_time are given, this workflow will automatically use --drmaa in its snakemake call.",
        required=False,
        type=int,
        default=None,
    )

    parser_workflow_cnv_calling.add_argument(
        "--drmaa_time",
        help="Maximum number of hh:mm per job. Suggested is 4:00.",
        required=False,
        type=str,
        default=None,
    )

    parser_workflow_cnv_calling.add_argument(
        "--drmaa_jobs",
        help="Number of jobs. Default is 50.",
        required=False,
        type=int,
        default=50,
    )
    parser_workflow_cnv_calling.set_defaults(func=workflow_cnv_calling.workflow_cnv_calling)

    # =========================================================================
    #  workflow reassignment
    # =========================================================================
    parser_workflow_reassignment = subparsers.add_parser(
        "workflow_reassignment",
        description=(
            "Complete CNV-calling snakemake-workflow to run on data from a single sequencing panel (bed-file)."
        ),
    )
    parser_workflow_reassignment.add_argument(
        "-w",
        "--workdir",
        help="Path to the snakemake workdir.",
        required=True,
        type=str,
    )
    parser_workflow_reassignment.add_argument(
        "-r",
        "--reference",
        help="Path to the genomic reference.",
        required=True,
        type=str,
    )
    parser_workflow_reassignment.add_argument(
        "-m",
        "--metafile",
        help="Path to the file containing the meta information. It is a .tsv of the scheme [panel] [bams.txt] [bed.bed]. \
            It aligns each desired panel name with the corresponding .bam files and the .bed file.",
        required=True,
        type=str,
    )

    parser_workflow_reassignment.add_argument(
        "-c",
        "--coverages",
        help="Output file path to coverages matrix (e.g. coverages.tsv).",
        required=True,
        type=str,
    )

    parser_workflow_reassignment.add_argument(
        "-b",
        "--bedfile",
        required=True,
        type=str,
        help="Output file path to BED file that will contain the union of all given BED files",
    )

    parser_workflow_reassignment.add_argument(
        "--cores",
        help="Number of used cores in snakemake workflow. Default is 32.",
        required=False,
        type=int,
        default=32,
    )

    parser_workflow_reassignment.add_argument(
        "--cluster_configfile",
        help="Path to the cluster config file.",
        required=False,
        type=str,
        default="",
    )

    parser_workflow_reassignment.add_argument(
        "--drmaa_mem",
        help="Number of megabytes used with drmaa. Suggested is 16000. If --drmaa_mem and --drmaa_time are given, this workflow will automatically use --drmaa in its snakemake call.",
        required=False,
        type=int,
        default=None,
    )

    parser_workflow_reassignment.add_argument(
        "--drmaa_time",
        help="Maximum number of hh:mm per job. Suggested is 4:00.",
        required=False,
        type=str,
        default=None,
    )
    parser_workflow_reassignment.add_argument(
        "--drmaa_jobs",
        help="Number of jobs. Default is 50.",
        required=False,
        type=int,
        default=50,
    )
    parser_workflow_reassignment.set_defaults(func=workflow_reassignment.workflow_reassignment)

    # =========================================================================
    #  visualize reassignment
    # =========================================================================
    parser_visualize_reassignmentment = subparsers.add_parser(
        "visualize_reassignment",
        description=("Interactive reassignment visualization"),
    )

    parser_visualize_reassignmentment.add_argument(
        "--host", default="0.0.0.0", help="Host to run server on, defaults to 0.0.0.0"
    )

    parser_visualize_reassignmentment.add_argument(
        "--port", type=int, default=8080, help="Port to run server on"
    )

    parser_visualize_reassignmentment.add_argument(
        "--debug", default=False, action="store_true", help="Whether or not to enable debugging"
    )

    parser_visualize_reassignmentment.add_argument(
        "--cache-dir", help="Optional path to cache directory, avoids repeating startup computation"
    )

    parser_visualize_reassignmentment.add_argument(
        "-m",
        "--metafile",
        help="Path to the file containing the meta information. It is a .tsv of the scheme [panel name] [path of bams.txt] [path of targets.bed]. \
            It aligns each desired panel name with the corresponding .bam files and the .bed file.",
        required=True,
        type=str,
    )

    parser_visualize_reassignmentment.add_argument(
        "-c",
        "--coverages",
        required=True,
        type=str,
        help="Matrix in tsv format containing coverage values. Is output of 'workflow_reassignment' step",
    )

    parser_visualize_reassignmentment.add_argument(
        "-b",
        "--bedfile",
        required=True,
        type=str,
        help="Path to BED file that will contain the union of all given BED files",
    )

    parser_visualize_reassignmentment.add_argument(
        "-d",
        "--new_panel_assignments_directory",
        help="Path to directory that will contain the output lists of paths to .bam-files, which were re-assigned to the given panels according to clustering.",
    )

    # parser_visualize_reassignmentment.add_argument(
    #    "-p",
    #    "--panels",
    #    required=True,
    #    type=str,
    #    help="Path to the directory that will contain the final lists of bam-files grouped by panels")

    # parser_visualize_reassignmentment.add_argument(
    #    "-a",
    #    "--batches",
    #    required=True,
    #    type=str,
    #    help="Path to the directory that will contain the final lists of bam-files grouped by batches",
    # )

    def viz_reassignment(args):
        """Helper function that launches the Dash webserver for reassignment visualization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            from .viz_reassignment import settings

            if args.cache_dir:
                settings.CACHE_DIR = args.cache_dir
            else:
                settings.CACHE_DIR = os.path.join(tmpdir, "cache")
            os.makedirs(settings.CACHE_DIR, exist_ok=True)
            settings.setup_paths(args)

            settings.BATCH_OUTPUT_PATH = args.new_panel_assignments_directory

            from .viz_reassignment.app import app  # noqa
            from .viz_reassignment import store

            if settings.CACHE_PRELOAD_DATA:
                store.load_all_data(settings.reassign_SETTINGS)
            from .viz_reassignment import ui_plots

            ## XXX
            us = settings.reassign_SETTINGS
            data = store.load_all_data(us)
            ui_plots.plot_clustermap_clustering_as_base64(
                data, *store.compute_acluster(us), store.compute_clustercoldict(us)
            )
            # ui_plots.plot_clustermap_batches_as_base64(
            #    data, store.compute_acluster(us), store.compute_clustercoldict(us)
            # )
            ## XXX
            app.run_server(host=args.host, port=args.port, debug=args.debug)

    parser_visualize_reassignmentment.set_defaults(func=viz_reassignment)

    return parser


def main():
    if not misc.execs_available(REQUIRED_EXECUTABLES):
        logger.error("Missing some executables. The program will eventually fail.")
    else:
        logger.debug("External executables present: %s", ", ".join(REQUIRED_EXECUTABLES))
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
