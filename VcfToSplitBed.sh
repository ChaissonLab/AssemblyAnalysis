#!/usr/bin/env bash
grep -v "^#" $1 | awk '{ a=length($4); b=length($5); if (a > b) print $1"\t"$2"\t"$2+a-b"\tdeletion\t"a-b"\t"$8;}' | bedtools sort > $2 & 
grep -v "^#" $1 | awk '{ a=length($4); b=length($5); if (b > a) print $1"\t"$2"\t"$2+b-a"\tinsertion\t"b-a"\t"$8;}' | bedtools sort > $3 &

wait
