#!/bin/bash
grna_pam=$1 ## eg:  AGAGCGTCCCATCTGGTAAGTCng
sample=$2 ## eg: BE4max_TCTT
awk ' {if($0~/^[0-9XY]+$/ )print $0}'  /path2refgnome/chr.txt | while read chr
do
    java -jar /path_calitas/calitas-1.0.jar SearchReference \
        -i $grna_pam \
        -I myguide  \
        -r /path2refgnome/GRCh38_${chr}.fa \
        -o /path_result/calitas/$sample/${chr}_myguide.hits.txt \
        --max-guide-diffs 3 \
        --max-pam-mismatches 0 \
        --max-gaps-between-guide-and-pam 2
done
