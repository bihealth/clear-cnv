PANEL=$1
BATCH=$2

#pip uninstall clearCNV -y
#python3 setup.py sdist bdist_wheel
#pip install -e .
rm /fast/work/users/vmay_m/workflow/FCC/CALLING/$PANEL\_$BATCH/results/cnv.calls
rm /fast/work/users/vmay_m/workflow/FCC/CALLING/$PANEL\_$BATCH/results/rscores.tsv
rm /fast/work/users/vmay_m/workflow/FCC/CALLING/$PANEL\_$BATCH/results/zscores.tsv
clearCNV workflow_cnv_calling \
-w /fast/users/vmay_m/work/workflow/FCC/CALLING/ \
-p $PANEL\_$BATCH \
-r /fast/work/projects/cubit/18.12/static_data/reference/GRCh37/hs37d5/hs37d5.fa \
-b /fast/work/groups/cubi/projects/2020-10-26_clear_cnv_analysis/data/untangling/run_001/batches/batch_$PANEL\_$BATCH.txt \
-d /fast/users/vmay_m/work/workflow/beds/nosex/$PANEL.bed \
-k /fast/work/projects/cubit/18.12/static_data/db/goldenpath/variable/GRCh37/wgEncodeCrgMapabilityAlign36mer.eq_1.bed \
--sensitivity 0.75 \
--cluster_configfile /fast/work/users/vmay_m/development/clearCNV/clearCNV/workflow/bih_cluster.json \
--drmaa_mem 32000 \
--drmaa_time 4:00 \
--cores 128

#-b /fast/users/vmay_m/work/workflow/bams/untangled/$PANEL\_bams.txt \
