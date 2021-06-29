import subprocess
import pathlib


def create_config(configpath, args):
    with open(configpath, "wt") as f:
        print("workdir: '" + args.workdir + "'", file=f)
        print("panelname: '" + args.panelname + "'", file=f)
        print("bedfile: '" + args.bedfile + "'", file=f)
        print("bamsfile: '" + args.bamsfile + "'", file=f)
        print("reference: '" + args.reference + "'", file=f)
        print("kmer_align: '" + args.kmer_align + "'", file=f)
        # options
        print("expected_artefacts: '" + str(args.expected_artefacts) + "'", file=f)
        print("sample_score_factor: '" + str(args.sample_score_factor) + "'", file=f)
        print("minimum_group_sizes: '" + str(args.minimum_group_sizes) + "'", file=f)
        print("zscale: '" + str(args.zscale) + "'", file=f)
        print("size: '" + str(args.size) + "'", file=f)
        print("del_cutoff '" + str(args.del_cutoff) + "'", file=f)
        print("dup_cutoff '" + str(args.dup_cutoff) + "'", file=f)
        print("trans_prob '" + str(args.trans_prob) + "'", file=f)
        print("plotregions '" + str(args.plot_regions) + "'", file=f)


def workflow_cnv_calling(args):
    configpath = (
        pathlib.Path(__file__).absolute().parent
        / pathlib.Path("workflow")
        / pathlib.Path("config_cnv_calling.yml")
    )
    create_config(configpath, args)

    arguments = [
        "snakemake",
        "--snakefile",
        str(
            pathlib.Path(__file__).absolute().parent
            / pathlib.Path("workflow")
            / pathlib.Path("Snakefile_cnv_calling")
        ),
        "--configfile",
        str(configpath),
        "--cores",
        str(args.cores),
    ]
    subprocess.check_call(arguments)

    # "panelname="+args.panelname,
    # "workdir="+args.workdir,
    # "bedfile="+args.bedfile,
    # "bamsfile="+args.bamsfile,
    # "reference="+args.reference,
    # "kmer_align="+args.kmer_align)
