#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import multiprocessing as mp
import math
from clearCNV import util

# def matchscores(panel, intsv_path, matchscores_path, expected_cnv_rate, fast=True):
def matchscores(args):
    # parse args
    intsv_path = args.coverages
    matchscores_path = args.matchscores
    expected_cnv_rate = args.expected_artefacts
    CORES = args.cores if args.cores else mp.cpu_count()
    fast = True if args.fast else False

    D0 = util.load_dataframe(intsv_path)
    D3a = util.normalize_within_exon(
        util.normalize_within_sample(
            util.filter_by_coverage(
                D0[[not (s.startswith("X") or s.startswith("Y")) for s in D0.index]]
            )
        )
    )
    # only take every i-th target into account to speed up
    if fast:
        X = D3a.T[D3a.index[:: min([len(D3a.index), math.ceil(len(D3a.index) / 1000)])]].T
    else:
        X = D3a

    preMatchscores = pd.DataFrame(np.zeros((len(X.columns), len(X.columns))))

    def collect_result(result):
        nonlocal preMatchscores
        preMatchscores.iloc[:, result[0]] = result[1]

    pool = mp.Pool(CORES)

    print(
        f"calc matchscores with {CORES} threads at an expected artefact rate "
        f"of {expected_cnv_rate}"
    )
    for i in range(len(X.columns)):
        pool.apply_async(util.matchscore, args=(i, X, expected_cnv_rate), callback=collect_result)
    pool.close()
    pool.join()

    Matchscores = preMatchscores.add(preMatchscores.T)
    Matchscores.index, Matchscores.columns = X.columns, X.columns
    Matchscores.to_csv(matchscores_path, sep="\t", header=True, index=True)
    return Matchscores
