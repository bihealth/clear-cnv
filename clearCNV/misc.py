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
    with open(args.outfile, "wt") as f:
        p1 = subprocess.Popen(
            ["bedtools", "merge", "-i", args.infile, "-c", "4", "-o", "collapse"],
            stdout=subprocess.PIPE,
        )
        p2 = subprocess.Popen(["sort", "-V", "-k1,1", "-k2,2"], stdin=p1.stdout, stdout=f)
        p1.stdout.close()
        p2.communicate()[0]


# =============================================================================
#  UNTANGLING - preparations


def merge_bedfiles(beds, bedfile):
    pathlib.Path(bedfile).absolute().parent.mkdir(parents=True, exist_ok=True)
    with open(bedfile, "wt") as f:
        p1 = subprocess.Popen(
            ["bedops", "--merge", *[line.rstrip("\n") for line in beds]], stdout=subprocess.PIPE,
        )
        p2 = subprocess.Popen(
            ["sort", "-V", "-k1,1", "-k2,2"], stdin=p1.stdout, stdout=subprocess.PIPE
        )
        p3 = subprocess.Popen(["bedtools", "merge"], stdin=p2.stdout, stdout=f)
        p3.communicate()[0]
        p1.stdout.close()
        p2.stdout.close()
    B = [l.rstrip("\n") for l in open(bedfile)]
    with open(bedfile, "wt") as g:
        for line in B:  # [::1+(int(len(B) / 1000))]:
            print(line, file=g)


def prepare_untangling(args):
    META = pd.read_csv(args.metafile, sep="\t")
    # prepare bamsfile
    pathlib.Path(args.bamsfile).absolute().parent.mkdir(parents=True, exist_ok=True)
    with open(args.bamsfile, "wt") as f:
        L = [[l.rstrip("\n") for l in open(bam)] for bam in META.iloc[:, 1]]
        for b in set([val for l in L for val in l]):
            print(b, file=f)
    merge_bedfiles(META.iloc[:, 2], args.bedfile)


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
            ["bedtools", "multicov", "-p", "-q", "20", "-bams", args.inbam, "-bed", args.bedfile],
            stdout=f,
        )


def merge_rtbeds(args):
    D = pd.read_csv(args.bedfile, sep="\t", header=None)
    dcols = D.shape[1]
    for i in range(len(args.rtbeds)):
        D[dcols + i] = pd.read_csv(args.rtbeds[i], sep="\t", header=None)[dcols]
    D.columns = ["chr", "start", "end", "gene"][:dcols] + [
        pathlib.Path(s).stem for s in args.rtbeds
    ]
    D.to_csv(args.coverages, sep="\t", index=False)
