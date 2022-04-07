# BE-dot
a comprehensive tool for sgRNA design and OT profile prediction
## Introduction:
BE-dot takes a given SNV as input, containing funtions about designing sgRNAs, predicting off-target profiles, and annotating off-target products. 
Besides providing precise correction schemes at DNA level,BE-dot can perform synonymous corrections at protein level by degeneracy. 
When predicting off-target profiles of single base editing systems, BE-dot integrates multiple tools such as Cas-OFFinder, CALITAS, CFD, CRISTA and uCRISPR. In addition, BE-dot can automatically analyze all possible editing products at off-target sites, and convert them into avinput format for functional annotation by ANNOVAR.
## Requirment
- python==3.6
- biopython==1.78
- pandas==1.1.5
- numpy==1.19.5

- Cas-OFFinder==3.0.0b3
- RNAstructure==6.4
## Usage:
### Overall direction
```
python3 BE-dot.py
usage: BE-dot.py

Program: BE-dot
              Version: 1.0
              Author: Zelu Wang
              Email: wangzl630@163.com


positional arguments:
  {designsgRNA_opt1,designsgRNA_opt2,OTprediction,OTannotation}
                        command
    designsgRNA_opt1    Design sgRNA for target SNV. #Option 1: input the SNV
                        and it's 100nt flanking sequence to design sgRNA
    designsgRNA_opt2    Design sgRNA for target SNV. #Option 2: input the
                        SNV's rsID to design sgRNA
    OTprediction        predict off-target profile
    OTannotation        Annotate off-target products

optional arguments:
  -h, --help            show this help message and exit
                 
```
### For designsgRNA_opt1:
```
usage: BE-dot.py designsgRNA_opt1 [-h] --jobID JOBID --upSeq <seq> --downSeq
                                  <seq> --mut MUT --wt WT
                                  [--codon_frame <int>] [-o OUTPUTPATH]

optional arguments:
  -h, --help            show this help message and exit
  --jobID JOBID         Give an ID to the sgRNA-design job
  --upSeq <seq>         50nt upstream to the SNV
  --downSeq <seq>       50nt downstream to the SNV
  --mut MUT             Mutation base
  --wt WT               Reference base
  --codon_frame <int>   Amino acid codon frame that SNV occupy, 1/2/3
  -o OUTPUTPATH, --outputPath OUTPUTPATH
                        Path of output file(s)

```
### For designsgRNA_opt2:
```
usage: BE-dot.py designsgRNA_opt2 [-h] --rsID RSID [-o OUTPUTPATH]

optional arguments:
  -h, --help            show this help message and exit
  --rsID RSID           SNV's ID in dbSNP
  -o OUTPUTPATH, --outputPath OUTPUTPATH
                        Path of output file(s)
```
### For OTprediction:
```
usage: BE-dot.py OTprediction [-h] -BE BE -grna GRNA -onpam ONPAM
                              [-genome <file>] [-mis MISMATCH_NUMBER]
                              [--DNAbulge DNABULGE] [--RNAbulge RNABULGE]
                              [-o OUTPUTPATH]

optional arguments:
  -h, --help            show this help message and exit
  -BE BE                BE to do OT-prediction
  -grna GRNA            Single gRNA sequence to analyse (20nt)
  -genome <file>, --target_genome <file>
                        Genome to search off-target sites
  -mis MISMATCH_NUMBER, --mismatch_number MISMATCH_NUMBER
                        Maximum mismatch number considered in aligning.
                        [default:3]
  --DNAbulge DNABULGE   Maximum DNA bulge number considered in aligning.
                        [default:1]
  --RNAbulge RNABULGE   Maximum RNA bulge number considered in aligning.
                        [default:1]
  -o OUTPUTPATH, --outputPath OUTPUTPATH
                        Path of output file(s)

```
### For OTannotation:
```
usage: BE-dot.py OTannotation [-h] [-BE BASE_EDITOR] [-i INPUT_FILE]
                              [-o OUTPUTPATH]

optional arguments:
  -h, --help            show this help message and exit
  -BE BASE_EDITOR, --base_editor BASE_EDITOR
                        BE to do OT-annotation
  -i INPUT_FILE, --input-file INPUT_FILE
                        Set input tsv file
  -o OUTPUTFILE, --outputFile OUTPUTFILE
                       

```
## Example
### designsgRNA_opt1:
```
python BE-dot.py designsgRNA_opt1 --jobID job001 --upSeq GCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAG
--downSeq GTCCCATCTGGTAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATG --mut C --wt T --codon_frame 1 --outputPath /usr/local/data/
```
### designsgRNA_opt2:
```
python BE-dot.py designsgRNA_opt2 --rsID rs80357410 --outputPath /path/data/
```
### OTprediction
```
python BE-dot.py OTprediction -BE BE-PLUS -grna AAATCTTAGAGCGTCCCATC -genome /usr/local/database/GRCh38.fa -o /usr/local/data/
```
### OTannotation
```
python BE-dot.py OTannotation -BE BE-PLUS -i /usr/local/data/BE-PLUS_AAAT.txt -o /usr/local/data/BE-PLUS_AAAT.avinput
```

