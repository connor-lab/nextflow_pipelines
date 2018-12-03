#!/usr/bin/env bash

prefix=$1

zcat ${prefix}.clean_1.fq.gz | awk '{if (NR%4 == 1) {print $1 "_" $2} else print}' | gzip -c > ${prefix}.renamed_1.fq.gz
zcat ${prefix}.clean_2.fq.gz | awk '{if (NR%4 == 1) {print $1 "_" $2} else print}' | gzip -c > ${prefix}.renamed_2.fq.gz
