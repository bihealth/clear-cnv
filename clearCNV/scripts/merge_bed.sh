bedtools merge -i $1 -c 4 -o collapse \
| sort -V -k1,1 -k2,2 \
> $2
