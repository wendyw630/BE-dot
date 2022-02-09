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
BE-dot.py [-h] [--jobID JOBID] [--upSeq <seq>] [--downSeq <seq>]
                 [--mutation MUTATION] [--wt WT] [--readingFrame <int>]
                 [--rsID RSID] [--CustomizeBE <file>] [--guide <seq>]
                 [--pam <seq>] [--target_genome <file>] [--BE_name BE_NAME]
                 [--OTs_file <file>]
optional arguments:
  -h, --help            show this help message and exit

### sgRNA design
 # Option 1: input the SNV and it's 50mer context sequence to design sgRNA:
  --jobID JOBID         give an ID to the sgRNA-design job
  --upSeq <seq>         25mer upstream to the SNV
  --downSeq <seq>       25mer downstream to the SNV
  --mutation MUTATION   the mutation base
  --wt WT               the ref. base
  --readingFrame <int>  the readingFrame that SNV occupy, 1/2/3

  # Option 2: input the SNV's rsID to design sgRNA:
  --rsID RSID

  # customized BE to add in BElist:
  --CustomizeBE <file>

### OT profile prediction
 # generate the cas-offinder input file:
  --guide <seq>         single gRNA sequence to analyse (20nt)
  --pam <seq>           pam to analyse
  --target_genome <file>
                        the genome to search off-target sites

### enumerate all the OT editing products, generating avinput file:
  --BE_name BE_NAME     the selected BE in former OT prediction step
  --OTs_file <file>     OTs file
                 
```




