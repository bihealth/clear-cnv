#!/usr/bin/env python
# coding: utf-8

import subprocess
import pandas as pd
import os
import pathlib
import tempfile

# =============================================================================
#  BED file merging


def merge_bedfile(args):
    f = open(args.outfile, "wt")
    p1 = subprocess.Popen(
        ["bedtools", "merge", "-i", args.infile, "-c", "4", "-o", "collapse"],
        stdout=subprocess.PIPE,
    )
    p2 = subprocess.Popen(["sort", "-V", "-k1,1", "-k2,2"], stdin=p1.stdout, stdout=f)
    p1.stdout.close()
    output = p2.communicate()[0]
    f.close()


# =============================================================================
#  ANNOTATIONS


def __gc_content(reference, bedfile, GC_annotation):
    # with path_out.open("wt") as f:
    with open(GC_annotation, "wt") as f:
        subprocess.check_call(["bedtools", "nuc", "-fi", reference, "-bed", bedfile], stdout=f)


def __mappability(bedfile, kmer_align, mappability_annotation):
    with open(mappability_annotation, "wt") as f:
        subprocess.check_call(["bedtools", "coverage", "-a", bedfile, "-b", kmer_align], stdout=f)


def __merge_annotations(mappability_annotation, GC_annotation, annotations):
    pd.read_csv(mappability_annotation, sep="\t", header=None).iloc[:, [0, 1, 2, 3, 6, 7]].join(
        pd.read_csv(GC_annotation, sep="\t").iloc[:, 5]
    ).to_csv(annotations, sep="\t", header=None, index=False)


def create_annotations_file(args):
    with tempfile.TemporaryDirectory() as tmp_dir:
        d = pathlib.Path(tmp_dir)
        GC_annotation = d / pathlib.Path("GC_annotation.bed")
        mappability_annotation = d / pathlib.Path("mappability_annotation.bed")
        __gc_content(args.reference, args.bedfile, GC_annotation)
        __mappability(args.bedfile, args.kmer_align, mappability_annotation)
        __merge_annotations(mappability_annotation, GC_annotation, args.annotations)


# =============================================================================
#  COVERAGES


def count_pair(args):
    with open(args.rtbed, "wt") as f:
        subprocess.call(
            ["bedtools", "multicov", "-p", "-q", "20", "-bams", args.inbam, "-bed", args.bedfile], stdout=f
        )


def merge_rtbeds(args):
    D = pd.read_csv(args.bedfile, sep="\t", header=None)
    for i in range(len(args.rtbeds)):
        D[5 + i] = pd.read_csv(args.rtbeds[i], sep="\t", header=None)[4]
    D.columns = ["chr", "start", "end", "gene"] + [
        os.path.basename(s).split(".")[0] for s in args.rtbeds
    ]
    D.to_csv(args.coverages, sep="\t", index=False)
