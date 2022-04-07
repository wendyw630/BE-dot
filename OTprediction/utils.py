#!/usr/bin/python
import os, sys, argparse
import pandas as pd
from pathlib import Path
import utils
from designsgRNA.baseEditorsTable import CBElist, ABElist
from OTprediction.CFD.CFD import calcCfdScore

# get POT list of gRNA based on Cas-OFFinder.
##### eg. :
# ref genome
# N*20NRG 1 1
# grnaNNN 3
## output:
def cas_offinder_input(be, grna, genome, mis, DNAbulge, RNAbulge):
    pam = None
    if be in CBElist:
        pam = CBElist[be][0]
    if be in ABElist:
        pam = ABElist[be][0]
    if pam:
        global pam
        sample_name = str(be) + grna[:4]
        global sample_name
        f = open("/root/project/data/cas_in_{0}.txt".format(sample_name), 'w')
        f.write(genome + '\n')
        f.write('NNNNNNNNNNNNNNNNNNNN' + str(pam).upper() + "\t" + str(DNAbulge) + "\t" + str(RNAbulge) + "\n")
        f.write(grna[:20] + 'N' * len(pam) + "\t" + str(mis))
        f.close()
        print("Search candidate off-target sites by Cas-OFFinder")
        os.system(
            "./cas-offinder /root/project/data/cas_in_{0}.txt C /root/project/data/cas_out_{0}.txt".format(sample_name))
    if not pam:
        exit("the input BE is not in BE-dot")


#######  get POT list of gRNA based on CALITAS
# eg:
#  AAAAAAAAAAACTTTTCCCGGnrg
### output:/root/project/data/calitas_${sample}.txt
def calitas_input(grna, mis, DNAbulge, RNAbulge):
    if pam:
        grna_pam = str(grna).upper() + str(pam).lower()
        bulge = DNAbulge + RNAbulge
        os.system("./calitas.sh {0} {1} {2} {3}".format(grna_pam, sample_name, mis, bulge))

    if not pam:
        exit("the input BE is not in BE-dot")


def OT_predict(be, grna):
    ###########   生成 BEdeepoff输入文件，并运行  ############
    if be in CBElist:
        BE_type = "CBE"
    if be in ABElist:
        BE_type = "ABE"
    infile = r'D:\00_project\predict\BE_dot_02\data\cas_out_{0}.txt'.format(sample_name)
    input_file = Path(infile)  ##eg: r'D:\00_project\predict\BE_dot_02\data\casoffinder_BE_PLUS_AAAT.txt'
    stem = input_file.stem
    output_dir = utils.safe_makedir(input_file.parent / 'eff_input')
    output_file = input_file.parent / 'eff_input' / (str(input_file.stem) + '_BEdeep.txt')

    if input_file:
        ##########  reform cas_outfile to BEdeepoff_infile  ##########
        df_inputs = pd.read_csv(input_file, sep='\t')
        df_inputs = df_inputs.drop('#Id', axis=1)
        cols = ['#Bulge_type', 'crRNA', 'DNA', 'Chromosome', 'Position', 'Direction',
                'Mismatches', 'Bulge_Size']
        df_inputs.columns = cols
        pam_len=len(BElist[be][0])
        crRNA = str(grna) + 'N' *pam_len
        DNA = str(grna) + df_inputs['DNA'][2][-(pam_len):]  ## just random pick the [2]row
        df_inputs = df_inputs.append(
            {'#Bulge_type': "X", 'crRNA': crRNA, 'DNA': DNA, 'Chromosome': 1, 'Position': 1, 'Direction': "+",
             'Mismatches': 0, 'Bulge_Size': 0}, ignore_index=True)
        df_inputs.to_csv(output_file, sep='\t', index=False)
        try:
            print('Predict with BEdeepoff')
            os.system("python ./bedeepoff.py -b {0} -i {1}".format(BE_type,output_file))
        except Exception as e:
            print(str(e))
        #################################################
        ###########  extract no-bulge OTs  ##########
        df_0bulge = df_inputs[df_inputs['#Bulge_type'] == "X"]

        ######   按行运行CFD   ######
        print('Predict with CFD')
        df_0bulge['CFD'] = df_0bulge.apply(lambda row: calcCfdScore(row['crRNA'], row['DNA']), axis=1)


        ###### uCRISPR score ######
        print('Predict with uCRISPR')
        os.system("./ucrispr.sh ")
        fucr = pd.read_csv('data/ucrispr.out', sep=' ', low_memory=False)
        df_0bulge= pd.concat([df_0bulge, fucr['uCRISPR']], axis=1)
        df_0bulge.to_csv(r'D:\00_project\predict\BE_dot_02\data\cas_out_{0}_0bulge.txt'.format(sample_name), sep='\t',
                         index=False)




