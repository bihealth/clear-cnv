import subprocess
import pathlib

from . import misc


def create_config(configpath, args):
    with open(configpath, "wt") as f:
        print("workdir: '" + str(args.workdir) + "'", file=f)
        print("metafile: '" + str(args.metafile) + "'", file=f)
        print("coverages: '" + str(args.coverages) + "'", file=f)
        print("cov_dir: '" + str(args.cov_dir) + "'", file=f)
        print("bedfile: '" + str(args.bedfile) + "'", file=f)
        print("bamsfile: '" + str(args.bamsfile) + "'", file=f)
        print("reference: '" + str(args.reference) + "'", file=f)

def workflow_reassignment(args):
    args.bamsfile = (pathlib.Path(args.workdir) / "allbams.txt").absolute()
    args.cov_dir = (pathlib.Path(args.workdir) / "covs").absolute()
    args.bedfile = pathlib.Path(args.bedfile).absolute()
    args.coverages = pathlib.Path(args.coverages).absolute()
    args.metafile = pathlib.Path(args.metafile).absolute()
    args.reference = pathlib.Path(args.reference).absolute()

    misc.prepare_reassignment(args)

    configpath = (
        pathlib.Path(__file__).absolute().parent
        / pathlib.Path("workflow")
        / pathlib.Path("config_reassignment.yml")
    )
    configpath.parent.mkdir(parents=True, exist_ok=True)
    create_config(configpath, args)
    arguments = [
        "snakemake",
        "--snakefile",
        str(
            pathlib.Path(__file__).absolute().parent
            / pathlib.Path("workflow")
            / pathlib.Path("Snakefile_reassignment")
        ),
        "--configfile",
        str(configpath),
        "--cores",
        str(args.cores),
    ]
    if args.cluster_configfile:
        arguments.append("--cluster-config")
        arguments.append(args.cluster_configfile)
    if args.drmaa_mem and args.drmaa_time:
        arguments.append("-p")
        arguments.append("--drmaa")
        arguments.append(
            str(
                f"\" --mem={args.drmaa_mem} --time={args.drmaa_time} --output={args.workdir}sge_log/%x-%j.log\""
            )
        )

    # How to create the dir the best way?
    subprocess.check_call(["mkdir", "-p", f"{args.workdir}sge_log"])
    subprocess.check_call(["mkdir", "-p", f"{args.workdir}covs"])

    subprocess.check_call(arguments)
