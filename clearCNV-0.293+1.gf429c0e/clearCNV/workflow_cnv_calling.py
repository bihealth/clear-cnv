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
        print("del_cutoff: '" + str(args.del_cutoff) + "'", file=f)
        print("dup_cutoff: '" + str(args.dup_cutoff) + "'", file=f)
        print("trans_prob: '" + str(args.trans_prob) + "'", file=f)
        print("plot_regions: '" + str(args.plot_regions) + "'", file=f)


def workflow_cnv_calling(args):
    configpath = (
        pathlib.Path(__file__).absolute().parent
        / pathlib.Path("workflow")
        / pathlib.Path("config_cnv_calling.yml")
    )
    configpath.parent.mkdir(parents=True, exist_ok=True)
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
    if args.cluster_configfile:
        arguments.append("--cluster-config")
        arguments.append(args.cluster_configfile)
    if args.drmaa_mem and args.drmaa_time and args.drmaa_jobs:
        arguments.append("-p")
        arguments.append("--drmaa")
        workdir_full_path = str(pathlib.Path(args.workdir).absolute())
        arguments.append(
            str(
                f'" --mem={args.drmaa_mem} --time={args.drmaa_time} --output={workdir_full_path}/sge_log/%x-%j.log"'
            )
        )
        arguments.append(f"--jobs {args.drmaa_jobs}")
    subprocess.check_call(arguments)

    # "panelname="+args.panelname,
    # "workdir="+args.workdir,
    # "bedfile="+args.bedfile,
    # "bamsfile="+args.bamsfile,
    # "reference="+args.reference,
    # "kmer_align="+args.kmer_align)
