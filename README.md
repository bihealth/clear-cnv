![Continuous Integration Status](https://github.com/bihealth/clear-CNV/workflows/CI/badge.svg)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/2eaafb57fbb74a46b918e9f58142c880)](https://www.codacy.com?utm_source=github.com&utm_medium=referral&utm_content=bihealth/clear-CNV&utm_campaign=Badge_Coverage)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/2eaafb57fbb74a46b918e9f58142c880)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bihealth/clear-CNV&amp;utm_campaign=Badge_Grade)
[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

# ClearCNV: Clinical sequencing panel CNV caller and visualizer

- Code Formatting: [black](https://github.com/psf/black)

## Installation

### conda

ClearCNV is available on conda: https://anaconda.org/bioconda/clearcnv

I'd recommend to create a conda env:

```mamba create -n clearcnv clearcnv -c conda-forge -c bioconda```

or

```conda create -n clearcnv clearcnv -c conda-forge -c bioconda```


Then clone this repo to your favorite location ```git clone git@github.com:bihealth/clear-cnv.git``` and `cd clear-cnv` into it. Lastly, activate the environment via ```mamba activate clearcnv``` or ```conda activate clearcnv```. Now you can run the commands listed below.

## Quickstart
#### Case one
You have one gene panel (e.g. called '/path/to/genepanel.bed') and a collection of aligned short reads (sample_xy0.bam, sample_xy1.bam, ..) and you want to call CNVs.
* Write a 'meta-file' like [this one](tests/testdata/meta.tsv).
  * Copy all full paths of your bam files to a txt-file e.g. '/path/to/bams.txt'.
  * Your '/path/to/meta.tsv' file would look like this: `genepanel\t/path/to/bams.txt\t/path/to/genepanel.bed`
* Use `clearCNV workflow_cnv_calling`. Type `clearCNV workflow_cnv_calling --help` to see how.
* Check the QC files to see if everything went well.
* Read more: https://github.com/bihealth/clear-cnv/blob/master/README.md#how-to-and-workflow

#### Case two
You have several panels and you're not really sure if the bam files are assigned correctly to each panel. You want the panels and batches separated and to call CNVs on each of them.
* Write a 'meta-file' like [this one](tests/testdata/meta.tsv).
  * Copy all full paths of your bam files that you think belong to panel 1 to a txt-file e.g. '/path/to/p1_bamfiles.txt'.
  * Copy all full paths of your bam files that you think belong to panel 2 to a txt-file e.g. '/path/to/p2_bamfiles.txt'. Do that for all panels. 
  * Your '/path/to/meta.tsv' file would look like this: `genepanel\t/path/to/bams.txt\t/path/to/genepanel.bed`
* Run `clearCNV workflow_reassignment`. Type `clearCNV workflow_reassignment --help` to see how.
* Run `clearCNV visualize_reassignment`. Type `clearCNV visualize_reassignment --help` to see how. You'll need to open the URL with your browser.
* After you ran each step in your browser, there will be a folder that contains all newly assigned batches. In each panel/batch you'll find a txt file that contains patchs to .bam files. These are your batches! Proceed with the `clearCNV workflow_cnv_calling` step for each batch. Type `clearCNV workflow_cnv_calling --help` to see how.
* Read more: https://github.com/bihealth/clear-cnv/blob/master/README.md#how-to-and-workflow

## Quick run checks and examples

### Sample reassignment:
#### Create all files

Execute the shell commamd (from within the cloned repo directory):
```clearCNV workflow_reassignment --workdir tests/testdata/ --reference tests/testdata/test_reassignment_ref.fa --metafile tests/testdata/test_reassign_meta.tsv --coverages tests/testdata/test_reassignment_coverages.tsv --bedfile tests/testdata/test_reassignment_union.bed --cores 2```

 - INPUT: working directory given by `--workdir`, the files given by `--reference` and `--metafile`.
 - OUTPUT: files created at `--coverages` and `--bedfile`. They are used in the next step.

If you want to create the necessary files for yourown data just edit the meta.tsv file analogously to the example at `clearCNV/tests/testdata/meta.tsv`, where you can add more rows for each targets file (BED-file). It is recommended to use absolute paths in the meta file.

Optionally, **drmaa** can be used, if the two flags are present:
`--drmaa_mem 1600 --drmaa_time 4:00`,
where drmaa is given 16 Gb memory per core and and four hours maximum running time.
Also, a cluster config file in .json format can be given with `--cluster_configfile config.json`

### Visualize sample reassignment:
#### Visualize and adjust the clusterings and final panel assignments

Execute the shell commamd (from within the cloned repo directory):
```clearCNV visualize_reassignment --metafile tests/testdata/meta.tsv --coverages tests/testdata/cov_reassignment.tsv --bedfile tests/testdata/reassignment_union.bed --new_panel_assignments_directory tests/testdata/panel_assignments```

 - INPUT: files given by `--metafile`, `--coverages` and `--bedfile`.
 - OUTPUT: files found in given directory `--new_panel_assignments_directory`.

### CNV calling

#### Match scores

At first, match scores are claculated. Go to the directory `clear-cnv/` and execute the shell command:

```clearCNV matchscores -p testpanel -c tests/testdata/cov.tsv -m tests/testdata/matchscores.tsv```

This creates a match score matrix which is used in the CNV calling step.

#### CNV calls

Now execute this shell command:

```clearCNV cnv_calling -p testpanel -c tests/testdata/cov.tsv -a tests/testdata/testpanel/analysis -m tests/testdata/matchscores.tsv -C tests/testdata/testpanel/results/cnv_calls.tsv -r tests/testdata/testpanel/results/rscores.tsv -z tests/testdata/testpanel/results/zscores.tsv -g 15 -u 3```

This creates the file `tests/testdata/testpanel/results/cnv_calls.tsv` which shows one called deletion. if you copy & paste this for your own data, please don't use the `-g 15 -u 3` configuration. We use these in here just to be able to work with a tiny example.

More files for analysis can now be found in `tests/testdata/testpanel/analysis`.

## HOW TO and WORKFLOW

clearCNV comprises of two major workflows and three major commads:

#### workflow

1) re-assignment (not necessary for CNV calling)
    
    a) `clearCNV workflow_reassignment`
    
    b) `clearCNV visualize_reassignment`
    
2) CNV calling
    
    a) `clearCNV workflow_cnv_calling`

#### preparations

Some files have to be acquired or created before these commands can be run:
1) re-assignment:
    
    a) For each sequencing panel a .bed file is needed following this [form](tests/testdata/panel1.bed).
    
    b) For each sequencing panel (or .bed-file containing all target informations) a simple list of the according .bam files is needed. An example can be found [here](tests/testdata/reassignment_p1_bamfiles.txt). Make sure to use absolute paths for this file.
    
    c) meta-file. This file is a tab-separated file and one example can be found [here](tests/testdata/meta.tsv). To avoid any confusion, we recommend using absolute paths here again.

2) CNV calling:

    a) A genome reference file. It must be the same that was used to create the read alignment files (.bam files).
    
    b) `workflow_cnv_calling` does CNV calling for each batch (or sequencing panel associated data set) separately. A text file with all .bam file paths for each batch and panel must be created. [Here](tests/testdata/test_reassignment_p1_bamfiles.txt) is an example showing only one .bam file path. Multiple paths are separated with a newline. This file is usually an output of `clearCNV visualize_reassignment`.
    
    c) The .bed-file for the sequencing panel for which this batch is put to CNV calling. An example can be found [here](tests/testdata/panel1.bed). Note that `gene` is optimally replaced with the real name of the exon, gene or target.
    
    d) A k-mer alignability file in .bed format. Such files can be downloaded from UCSC (e.g. for Hg19 [here](http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability)). A k-mer mappability track can also be created for example using [GenMap](https://github.com/cpockrandt/genmap). In both cases the resulting Wig or BigWig files need to be converted to .bed to be used by clearCNV.

#### NOTES
The chromosome name scheme in the reference and .bed-file should be of the forms: ChrX, chrX, X or Chr1, chr1, 1.

CNV calling on chr X or chr Y: clearCNV automatically determines the copy number of the gonosomes. If your panel targets only a single gene per chromosome, then it is better to delete according targets from the original .bed file to exclude them. It is necessary to have about double as many samples in your data set to enable meaningful CNV calling on the X or Y chromosomes with roughly equally many women and men in the samples.

If you do sample re-assignment on your own data, followed by CNV-calling, then only one metafile, one coverages file, and one bedfile will be used. This means that `--metafile`, `--coverages` and `--bedfile` are given the same file paths in both workflow steps `clearCNV workflow_reassignment` and `clearCNV visualize_reassignment` of clearCNV. The coverages file can not be re-used for the CNV calling steps.

### Running Checks

Checks are automatically run on the `master` branch and pull requests.
Unit and integration tests are based on pytest and formatting is enforced with black.


```bash
$ make test
```
