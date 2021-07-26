#!/usr/bin/env python

import contextlib
import pathlib
import subprocess
import tempfile
import typing

from logzero import logger
import pandas as pd


def execs_available(execs: typing.List[str]) -> None:
    ok = True
    for exec in execs:
        res = subprocess.run(["which", exec], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if res.returncode != 0:
            logger.warn("Could not find executable %s via `which %s`", exec, exec)
            ok = False
    return ok


def merge_bedfile(args):
    # TODO: used anywhere besides tests?
    # Yes! Each time target counts are created from bam files!
    cmd_merge = ["bedtools", "merge", "-i", args.infile, "-c", "4", "-o", "collapse"]
    cmd_sort = ["sort", "-V", "-k1,1", "-k2,2"]
    with open(args.outfile, "wt") as f:
        with subprocess.Popen(cmd_merge, stdout=subprocess.PIPE, shell=False) as proc_merge:
            with subprocess.Popen(
                cmd_sort, stdin=proc_merge.stdout, stdout=f, shell=False
            ) as proc_sort:
                proc_merge.wait()
                proc_sort.wait()


def _popen(cmd, stdin=None, stdout=None):
    """Simple wrapper for Popen to shorten things and in particular disables using the shell."""
    stdout = stdout or subprocess.PIPE
    return subprocess.Popen(cmd, stdin=stdin, stdout=stdout, shell=False)


# =============================================================================
#  reassignment - preparations


def merge_bedfiles(beds, bedfile):
    pathlib.Path().absolute().parent.mkdir(parents=True, exist_ok=True)

    cmd_merge = ["bedops", "--merge", *map(str.rstrip, beds)]
    cmd_sort = ["sort", "-V", "-k1,1", "-k2,2"]
    cmd_merge2 = ["bedtools", "merge"]

    with contextlib.ExitStack() as stack:
        # Open output file and start process pipeline.
        f = stack.enter_context(open(bedfile, "wt"))
        p_merge = _popen(cmd_merge)
        p_sort = _popen(cmd_sort, p_merge.stdout)
        p_merge2 = _popen(cmd_merge2, p_sort.stdout, f)
        # Wait for completion.
        p_merge.wait()
        p_sort.wait()
        p_merge2.wait()


def prepare_reassignment(args):
    META = pd.read_csv(args.metafile, sep="\t")
    # prepare bamsfile
    pathlib.Path(args.bamsfile).absolute().parent.mkdir(parents=True, exist_ok=True)
    f = open(args.bamsfile, "wt")
    L = [[l.rstrip("\n") for l in open(bam)] for bam in META.iloc[:, 1]]
    for b in set([val for l in L for val in l]):
        print(b, file=f)
    f.close()
    merge_bedfiles(list(META["bedfiles"]), args.bedfile)


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
