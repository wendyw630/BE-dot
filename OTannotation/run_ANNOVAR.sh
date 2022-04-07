#!/bin/bash
sample=$1
/usr/local/software/annovar/annotate_variation.pl \
 -out /usr/local/data/${sample} \
 -build hg38 \
 /usr/local/data/${sample}.avinput \
 /usr/local/software/annovar/humandb/
