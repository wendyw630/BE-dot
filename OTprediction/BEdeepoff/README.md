### Statement
all the files in this folder are referenced from 
"Zhang C, Wang D, Qi T, et al. BEdeepoff: an in silico tool for off-target prediction of ABE and CBE base editors. bioRxiv, 2021"

### System requirements
The code were tesed on Linux and Mac OS systems.

The required software/packages are:
* python>=3.6.5
* numpy=1.19.1
* Pytorch=1.6.0
* scikit-learn
* biopython=1.78
* pandas


### Help
- bedeepon.py On-target models for ABEmax and AncBE4max

  -b {ABE,CBE} set base editor model
  
  -i INPUT_FILE set input tsv file

- bedeepoff.py Off-target models for ABEmax and AncBE4max

  -b {ABE,CBE} Set base editor model

  -i INPUT_FILE Set input tsv file

### Example Call
- On-target ABE prediction:
```bash
python bedeepon.py -b ABE -i ./Example/On-target.tsv
```
- Off-target ABE prediction for Cas-OFFinder output:
```bash
python bedeepoff.py -b ABE -i ./cas-offinder/Cas-Offinder_Output_Example.txt
```
- Off-target ABE prediction for CRISPRitz output:
```bash
python bedeepoff.py -b ABE -i ./cas-offinder/CRISPRitz_Output_Example.txt
```
