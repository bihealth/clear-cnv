import subprocess
import os
import pathlib

from . import cnv_arithmetics as ca
from . import matchscores
from . import cnv_calling
from . import visualize_scores


def create_config(configpath, args):
    with open(configpath, "wt") as f:
        print("workdir: '" + args.workdir + "'", file=f)
        print("panelname: '" + args.panelname + "'", file=f)
        print("bedfile: '" + args.bedfile + "'", file=f)
        print("bamsfile: '" + args.bamsfile + "'", file=f)
        print("reference: '" + args.reference + "'", file=f)
        print("kmer_align: '" + args.kmer_align + "'", file=f)


def workflow(args):
    configpath = (
        os.path.dirname(os.path.abspath(__file__))
        / pathlib.Path("workflow")
        / pathlib.Path("config.yml")
    )
    create_config(configpath, args)
    subprocess.check_call(
        [
            "snakemake",
            "--snakefile",
            str(
                os.path.dirname(os.path.abspath(__file__))
                / pathlib.Path("workflow")
                / pathlib.Path("Snakefile")
            ),
            "--configfile",
            str(configpath),
        ]
    )

    # "panelname="+args.panelname,
    # "workdir="+args.workdir,
    # "bedfile="+args.bedfile,
    # "bamsfile="+args.bamsfile,
    # "reference="+args.reference,
    # "kmer_align="+args.kmer_align)
