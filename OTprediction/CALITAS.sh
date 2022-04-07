#!/bin/bash
grna_pam=$1 ## eg:  AGAGCGTCCCATCTGGTAAGTCng
sample=$2 ## eg: BE4max_TCTT
max_mis=$3 ## eg: 3
max_bulge=$4 ## eg: 2
awk ' {if($0~/^[0-9XY]+$/ )print $0}'  /usr/local/database/GRCh38/ref_split/chr.txt | while read chr
do
    java -jar /usr/local/tool/calitas/calitas-1.0.jar SearchReference \
        -i $grna_pam \
        -I myguide  \
        -r /usr/local/database/GRCh38/ref_split/GRCh38_${chr}.fa \
        -o /usr/local/predict/result/calitas/$sample/${chr}_myguide.hits.txt \
        --max-guide-diffs $max_mis \
        --max-pam-mismatches 0 \
        --max-gaps-between-guide-and-pam $max_bulge
done

######  writting   #######
cd  /root/project/data/
cat /dev/null > calitas_$sample.txt
echo -e '##padded_guide\tpadded_target\tchromosome\tcoordinate_start\tstrand\tscore' > calitas_${sample}.txt

for file in  /usr/local/predict/result/calitas/$sample/*txt
do
    file_l=$(cat $file | wc -l )
    if [ $file_l -gt 1 ]
    then
        awk 'BEGIN{FS="\t";OFS="\t"}NR>1{print $22,$24,$4,$5,$7,$16}' $file >> calitas_${sample}.txt
    fi
done

awk 'BEGIN{FS="\t"; OFS=" "}NR>1{if($1~/^\w+$/  && $2~/^\w+$/) print $1,$2}' calitas_${sample}.txt > calitas_${sample}.txt
