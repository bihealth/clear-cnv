#!/usr/bin/env python
# coding: utf-8

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description="Merges coverages per sample (.rtbed files) into one dataframe.")
parser.add_argument("-b", "--bed", help="bed-file of sequencing panel",required=True, type=str)
parser.add_argument("-r", "--rtbeds", help="Paths to all .rtbed files.", required=True, nargs="+", type=str)
parser.add_argument("-c", "--coverages", help="Output file of concatenated coverages", required=True, type=str)
args = parser.parse_args()

exome_path        = args.bed
rtbed_paths       = args.rtbeds
final_file        = args.coverages

def merge_rtbeds(bedfile, rtbed_paths, coverages):
    D = pd.read_csv(bedfile, sep='\t',header=None)
    for i in range(len(rtbed_paths)):
        D[5+i] = pd.read_csv(rtbed_paths[i], sep='\t', header=None)[4]
    D.columns=['chr','start','end','gene']+[os.path.basename(s).split('.')[0] for s in rtbed_paths]
    D.to_csv(coverages, sep='\t',index=False)
