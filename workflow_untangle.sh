#pip uninstall clearCNV -y
#python3 setup.py sdist bdist_wheel
#pip install -e .

#clearCNV workflow_untangle \
#-w /fast/work/users/vmay_m/workflow/FCC/UNTANGLING/run_002 \
#-r /fast/work/projects/cubit/18.12/static_data/reference/GRCh37/hs37d5/hs37d5.fa \
#-m /fast/work/users/vmay_m/workflow/FCC/UNTANGLING/run_001/meta.tsv
#--cores 16

clearCNV workflow_untangle \
-w /fast/work/users/vmay_m/workflow/FCC/UNTANGLING/run_001 \
-r /fast/work/projects/cubit/18.12/static_data/reference/GRCh37/hs37d5/hs37d5.fa \
-m /fast/work/users/vmay_m/workflow/FCC/UNTANGLING/run_001/meta.tsv \
--cluster_configfile /fast/work/users/vmay_m/development/clearCNV/clearCNV/workflow/bih_cluster.json \
--drmaa_mem 32000 \
--drmaa_time 4:00 \
--cores 256
