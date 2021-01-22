PANEL=$1

#pip uninstall clearCNV -y
#python3 setup.py sdist bdist_wheel
#pip install -e .

#rm /fast/work/users/vmay_m/workflow/FCC/CALLING/$PANEL/results/cnv.calls
#rm /fast/work/users/vmay_m/workflow/FCC/CALLING/$PANEL/results/rscores.tsv
#rm /fast/work/users/vmay_m/workflow/FCC/CALLING/$PANEL/results/zscores.tsv

rm /fast/work/users/vmay_m/workflow/FCC/CALLING/$PANEL/analysis/z_scores_extended_*.html
rm /fast/work/users/vmay_m/workflow/FCC/CALLING/$PANEL/analysis/ratio_scores_extended_*.html

clearCNV workflow_cnv_calling \
-w /fast/users/vmay_m/work/workflow/FCC/CALLING/ \
-p $PANEL \
-r /fast/work/projects/cubit/18.12/static_data/reference/GRCh37/hs37d5/hs37d5.fa \
-b /fast/work/groups/cubi/projects/2020-10-26_clear_cnv_analysis/data/untangling/run_001/$PANEL\_clustered_bams.txt \
-d /fast/users/vmay_m/work/workflow/beds/nosex/$PANEL.bed \
-k /fast/work/projects/cubit/18.12/static_data/db/goldenpath/variable/GRCh37/wgEncodeCrgMapabilityAlign36mer.eq_1.bed \
--sensitivity 0.75 \
--minimum_group_sizes 25 \
--minimum_sample_score 0.15 \
--cluster_configfile /fast/work/users/vmay_m/development/clearCNV/clearCNV/workflow/bih_cluster.json \
--drmaa_mem 64000 \
--drmaa_time 4:00 \
--cores 128

#-b /fast/users/vmay_m/work/workflow/bams/untangled/$PANEL\_bams.txt \
