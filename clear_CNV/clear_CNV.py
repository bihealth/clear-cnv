#!/usr/bin/env python
# coding: utf-8

# main parsing
def get_parser():
    import argparse
    import util
    import cnv_arithmetics as ca
    import matchscores
    import cnv_calling
    import visualize_scores
    parser = argparse.ArgumentParser(prog='clear_CNV',  usage='%(prog)s [options]', description='clear-CNV can compute matchscores, CNV-calls and visualizations.')
    subparsers = parser.add_subparsers()

    parser_matchscores = subparsers.add_parser('matchscores', description="Matchscore calculation script.")
    parser_matchscores.add_argument("-p", "--panel",              help="Name of the data set (or panel)",required=True, type=str)
    parser_matchscores.add_argument("-c", "--coverages",          help="Coverages file in tsv format",required=True, type=str)
    parser_matchscores.add_argument("-m", "--matchscores",        help="Output matchscores.tsv file", required=True, type=str)
    parser_matchscores.add_argument("-x", "--expected_artifacts", help="Expected ratio of CNVs or artifacs in target fragment counts", required=False, type=float, default=0.02)
    parser_matchscores.add_argument("--cores",                    help="Number of cpu cores used in parallel processing. Default: determined automatically.", required=False, type=int,   default=0)
    parser_matchscores.set_defaults(func=matchscores.matchscores)

    parser_cnv_calling = subparsers.add_parser('cnv_calling',description="CNV calling script. Output is a single file in tsv format containing a list of CNV calls sorted by score. Some quality control plots are added to the analysis directory in the process.")
    parser_cnv_calling.add_argument("-p", "--panel",               help="Name of the data set(or panel)",                                                      required=True,  type=str)
    parser_cnv_calling.add_argument("-c", "--coverages",           help="Coverages file in tsv format",                                                        required=True,  type=str)
    parser_cnv_calling.add_argument("-a", "--analysis_directory",  help="Path to the directory, where analysis files are stored",                              required=True,  type=str)
    parser_cnv_calling.add_argument("-m", "--matchscores",         help="matchscores.tsv file generated with matchscores.py",                                  required=True,  type=str)
    parser_cnv_calling.add_argument("-C", "--cnv_calls",           help="Output cnv.calls file formatted in tsv format",                                       required=True,  type=str)
    parser_cnv_calling.add_argument("-r", "--ratio_scores",        help="Output ratio scores file in tsv format. Best kept together with cnv.calls",           required=True,  type=str)
    parser_cnv_calling.add_argument("-z", "--z_scores",            help="Output z-scores file in tsv format. Best kept together with cnv.calls",               required=True,  type=str)
    parser_cnv_calling.add_argument("-x", "--expected_artifacts",  help="Expected ratio of CNVs or artifacs in target fragment counts",                        required=False, type=float, default=0.02)
    parser_cnv_calling.add_argument("-u", "--minimum_sample_score",help="A lower threshold results in better fitting, but smaller calling groups",             required=False, type=float, default=0.15)
    parser_cnv_calling.add_argument("-g", "--minimum_group_sizes", help="Minimum group size per CNV calling group per match scores",                           required=False, type=int,   default=33)
    parser_cnv_calling.add_argument("-s", "--sensitivity",         help="A higher sensitivity results in more CNV calls. Can only be 0.0 <= sens <= 1.0",      required=False, type=float, default=0.7)
    parser_cnv_calling.add_argument("--cores",                     help="Number of cpu cores used in parallel processing. Default: determined automatically.", required=False, type=int,   default=0)
    parser_cnv_calling.set_defaults(func=cnv_calling.cnv_calling)

    parser_visualize = subparsers.add_parser('visualize',description="The visualization script creates html files containing heatmap-like matrices containing the ratios aligned with mappability, GC-content and target size so that CNVs can be visually identified and evaluated easily.")
    parser_visualize.add_argument("-a", "--analysis_directory", help="Path to the directory, where analysis files are stored", required=True, type=str)
    parser_visualize.add_argument("-r", "--ratio_scores",       help="Ratio scores file in tsv format, generated in cnv_calling.py", required=True, type=str)
    parser_visualize.add_argument("-z", "--z_scores",           help="Z-scores file in tsv format, generated in cnv_calling.py", required=True, type=str)
    parser_visualize.add_argument("-n", "--annotated",          help="CALL_GROUPS.py analysis_directory z_scores.tsv ratio_scores.tsv annotations.bed", required=True, type=str)
    parser_visualize.add_argument("-s", "--size",               help="Rough number of targets in each visualization.", required=False, type=int, default=1000)
    parser_visualize.set_defaults(func=visualize_scores.visualize)
    return parser

def main():
    from datetime import datetime
    print("Time at start:", datetime.now())
    import argparse
    args = get_parser().parse_args()
    args.func(args)
    print("Terminated at:", datetime.now())

if __name__ == '__main__':
    main()
