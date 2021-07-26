![Continuous Integration Status](https://github.com/bihealth/clear-CNV/workflows/CI/badge.svg)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/2eaafb57fbb74a46b918e9f58142c880)](https://www.codacy.com?utm_source=github.com&utm_medium=referral&utm_content=bihealth/clear-CNV&utm_campaign=Badge_Coverage)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/2eaafb57fbb74a46b918e9f58142c880)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bihealth/clear-CNV&amp;utm_campaign=Badge_Grade)
[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

# ClearCNV: Clinical sequencing panel CNV caller and visualizer

- Code Formatting: [black](https://github.com/psf/black)

## Developer Documentation

### HOW TO - panel reassignment
#### 1. Create all files to do reassignment

Go to the directory `clear-cnv/` and execute the shell commad:
```clearCNV workflow_reassign --workdir tests/testdata/ --reference tests/testdata/test_reassignment_ref.fa --metafile tests/testdata/test_reassign_meta.tsv --coverages tests/testdata/test_reassignment_coverages.tsv --bedfile tests/testdata/test_reassignment_union.bed --cores 2```

INPUT: working directory given by `--workdir`, the files given by `--reference` and `--metafile`.
OUTPUT: files created at `--coverages` and `--bedfile`. They are used in the next step.

If you want to create the necessary files for your own data just edit the meta.tsv file analogously to the example at `clearCNV/tests/testdata/meta.tsv`, where you can add more rows for each targets file (BED-file). It is recommended to use absolute paths in the meta file.

Optionally, **drmaa** can be used, if the two flags are present:
`--drmaa_mem 1600 --drmaa_time 4:00`,
where drmaa is given 16 Gb memory per core and and four hours maximum running time.
Also, a cluster config file in .json format can be given with `--cluster_configfile config.json`

#### 2. Visualize and adjust the clusterings and final panel assignments

Run the following shell command from `clear-cnv/`:
```clearCNV visualize_reassign --metafile tests/testdata/meta.tsv --coverages tests/testdata/cov_reassignment.tsv --bedfile tests/testdata/reassignment_union.bed --new_panel_assignments_directory tests/testdata/panel_assignments```

INPUT: files given by `--metafile`, `--coverages` and `--bedfile`.
OUTPUT: files found in given directory `--new_panel_assignments_directory`.

#### NOTE
If you do sample re-asignment on your own data, followed by CNV-calling, then only one metafile, one coverages file and one bedfile will be used. This means that `--metafile`, `--coverages` and `--bedfile` are given the same file paths in both workflow steps `clearCNV workflow_reassign` and `clearCNV visualize_reassign` of clearCNV.

### Running Checks

Checks are automatically run on the `master` branch and pull requests.
Unit and integration tests are based on pytest and formatting is enforced with black.


```bash
$ make test
```
