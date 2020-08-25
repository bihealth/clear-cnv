common_cnvs=$2
calls=$1
annotated=$3

gunzip $common_cnvs

echo bedtools intersect -wb -a <(awk '{print $2,$3,$4,$6,$7,$1}' $common_cnvs \
  | sed 's/ /\t/g') -b $calls > $annotated

echo "cchr  cstart  cend  cabb  cref	caccession	chr start end	abb gene  ratio zscore  sample  sampleScore" > $annotated

bedtools intersect -wb -a <(awk '{print $2,$3,$4,$6,$7,$1}' $common_cnvs \
 | sed 's/ /\t/g') -b $calls >> $annotated

gzip $common_cnvs

#bedtools intersect -v -a <(awk '{print $2,$3,$4,$6}' $common_cnvs | sed 's/ /\t/g') -b $calls >> $annotated
