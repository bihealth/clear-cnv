paste \
<(awk '{{print $1,$2,$3,$4,$7,$8}}' $1 | sed 's/ /\t/g') \
<(awk '{{print $6}}' $2 | sed 's/ /\t/g' | tail -n +2) \
> $3
