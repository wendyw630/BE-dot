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

```

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




