#!/usr/bin/python
######### author: wangzelu ###########

#####  sgRNA==20nt

from designsgRNA.baseEditorsTable import CBElist,  ABElist
import os
import sys
import itertools
import re
import pandas as pd
from Bio.Seq import Seq
from pathlib import Path

def OT_anno(be,infile,outfile):
    if be in CBElist:
        BElist=CBElist
    elif be in ABElist:
        BElist=ABElist
    try:
        windowl = BElist[be][1]
        windowu = BElist[be][2]
        infile=pd.read_table(infile,sep='\t')
        fout = open(outfile, 'w')
        for index,row in infile.iterrows():

            OT_seq = row['DNA'].upper()
            pat = '-'
            OT_seq = re.sub(pat, '', OT_seq)  ### for row['#bulge_type']=='RNA'

            C_rpos_lt = []
            ###divide conditions by 'bulge_type'
            #if row['#Bulge type']=='DNA' or row['#Bulge type']=='X':

            genome_windowl = windowl
            genome_windowu = windowu

            for pos in range(genome_windowl - 1, genome_windowu):
                if OT_seq[pos] == 'C':### rpos: 0-based
                    C_rpos_lt.append(pos)

            C_rpos_comb_lt = []
            if len(C_rpos_lt) > 0:
                for n in range(1, len(C_rpos_lt) + 1):
                    for C_rpos_comb in itertools.combinations(C_rpos_lt, n):
                        C_rpos_comb_lt.append(C_rpos_comb)


            if row['Direction']=="+":
                ## output
                for C_rpos_comb in C_rpos_comb_lt:
                    inpos1=(row['Location']+1)+min(C_rpos_comb) ### +1:in cas-offinder result, 'Position' should +1
                    inpos2=(row['Location']+1)+max(C_rpos_comb)
                    ### C_apos_comb
                    C_apos_comb = []
                    for rpos in C_rpos_comb:
                        C_apos_comb.append(rpos + 1)### apos=rpos+1

                    ### ref_frag && alt_frag
                    ref_frag=OT_seq[min(C_rpos_comb):max(C_rpos_comb)+1]

                    alt_frag_lt=[i for i in OT_seq[0:]]
                    for rpos in C_rpos_comb:
                        alt_frag_lt[rpos]='T'
                    alt_frag=''.join(alt_frag_lt[min(C_rpos_comb):max(C_rpos_comb)+1])
                    m = re.match(r'^([0-9XY]+)', row['Chromosome'])  ### wzl_change2
                    outchromosome = m.group(1)  ### wzl_change2
                    fout.write(outchromosome+"\t"+str(inpos1)+"\t"+str(inpos2)+"\t"+
                               ref_frag+"\t"+alt_frag+"\t"+str(OT_seq)+","+str(C_apos_comb)+"\n")


            if row['Direction']=='-':
                ## reverse_complement
                rev_comp_OT=str(Seq(OT_seq).reverse_complement())
                ## output
                for C_rpos_comb in C_rpos_comb_lt:

                    inpos1 = (row['Location']+1)+(len(OT_seq)-1) - max(C_rpos_comb) ### len(OT_seq): considering len(row['DNA'])>23
                    ## pre: row['Position']+len(OT_seq) - max(C_rpos_comb)+1
                    inpos2 = (row['Location']+1)+(len(OT_seq)-1) - min(C_rpos_comb)
                    ### C_apos_comb
                    C_apos_comb = []
                    for rpos in C_rpos_comb:
                        C_apos_comb.append(rpos + 1)

                    ### ref_frag && alt_frag
                    ref_frag = rev_comp_OT[len(OT_seq)-1-max(C_rpos_comb): len(OT_seq)-1-min(C_rpos_comb)+1]
                    alt_frag_lt=[i for i in rev_comp_OT[0:]]
                    for apos in C_apos_comb:
                        alt_frag_lt[-apos]='A'  ## rpos+1==apos
                    alt_frag = ''.join(alt_frag_lt[len(OT_seq)-1-max(C_rpos_comb): len(OT_seq)-1-min(C_rpos_comb)+1])

                    m = re.match(r'^([0-9XY]+)', row['Chromosome'])  ### wzl_change2
                    outchromosome = m.group(1)  ### wzl_change2
                    fout.write(outchromosome + "\t" + str(inpos1) + "\t" + str(inpos2) + "\t" +
                               ref_frag + "\t" + alt_frag + "\t" + rev_comp_OT+ "," + "minus"+str(C_apos_comb) + "\n")
            #if row['#Bulge type']=='RNA':
        fout.close()
        ###############################################
        ##############  run ANNOVAR  ##############
        anno_infile = Path(outfile)
        sample_name=str(anno_infile.stem)
        os.system("./OTannotation/run_ANNOVAR.sh {0}".format(sample_name))
    except Exception as e:
        print(str(e))


'''
if __name__=='__main__':
    OT_anno('BE-PLUS',r'/ANNOVAR/casoffinder_4Xots_BE_PLUS_AAAT.txt',
         r'/ANNOVAR/casoffinder_4Xots_BE_PLUS_AAAT.avinput')
    print("done")
'''
