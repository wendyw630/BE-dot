from designsgRNA.baseEditorsTable import CBElist,  ABElist, BEletter
from designsgRNA.snpClass import snp_define,supo_snp_define
from designsgRNA.fetchrsID import get_rsIDsnp


from Bio.Seq import Seq
import csv
import re
import copy
###############################################################

hive_dic={
    "BE4(max)":"BE4",
    "eA3A-BE3":"eA3A",
    "CP-CBEmax variants":"BE4-CP1028",
    "evoAPOBEC1-BE4max":"evoAPOBEC",
    "AID-BE4max":"AID",
    "CDA1-BE4max":"CDA",
    "ABE7.10":"ABE",
    "CP-ABEmax variants":"ABE-CP1040"
}

def hive_score(be, seq5,seq3,mut, grna): #
    pred_d=None
    if be in CBElist:
        BElist=CBElist
    else:
        BElist=ABElist
    if be in hive_dic:
        hive_be=hive_dic[be]
        full_seq=seq5+mut+seq3
        grna_seg=grna[-(BElist[be][1]-1):]

        if re.search(grna_seg,full_seq): ## rc_flag=False
            pass
        else:  ## rc_flag=True
            full_seq = Seq(full_seq).reverse_complement()
        try:
            seg_end_cor = re.search(grna_seg, full_seq).span()[1]
            seq=full_seq[seg_end_cor-40 : seg_end_cor+10]
            print("seq to writein hive",seq)
        except Exception as e:
            print(str(e))

        #be_efficiency_model.init_model(base_editor=hive_be, celltype='mES')
        #pred_d = str(be_efficiency_model.predict(seq))
        pred_d="k"
    else:
        pred_d="unknown"

    return pred_d

def cut_3nlen(snp):  ### from dot_BEsinglesyn3.py
    f1 = lambda x, i: x%((x // 3 - 1) * 3 + (i-1))## x: seq_len; i: readingframe
    f2 = lambda x,i: (x // 3 - 1) * 3 + (3-i)
    if snp.readingframe==1:
        snp.seq5=snp.seq5[f1(25,1):]
        snp.seq3=snp.seq3[:f2(25,1)]
    if snp.readingframe==2:
        snp.seq5=snp.seq5[f1(25,2):]
        snp.seq3=snp.seq3[:f2(25,2)]
    if snp.readingframe==3:
        snp.seq5=snp.seq5[f1(25,3):]
        snp.seq3=snp.seq3[:f2(25,3)]
    return snp
def synonymous12(snp,rc_flag,under_editing_cor,BElist,BE):
    ## snp.wt,mut,readingFrame,seq5,seq3
    synonymous12_flag=False
    #### case 1, case2
    ## TRANSLATE the whole seq with ref_'T' in under_editing_seq
    ## TRANSLATE the whole seq with editresult_'T' in under_editing_seq
    ### under_editing_seq, after_editing_seq:
    if under_editing_cor[0] !=0:
        under_editing_seq=snp.seq5[under_editing_cor[0]:]+snp.mut+snp.seq3[:under_editing_cor[1]]
    else:
        under_editing_seq=snp.mut+snp.seq3[:under_editing_cor[1]]
    after_editing_seq=None;after_full_seq=None 
    if BElist[BE][3]=='C':
        after_editing_seq=under_editing_seq.replace('C','T')
        #print("after_editing_seq:\t", after_editing_seq)
    if BElist[BE][3]=='A':
        after_editing_seq=under_editing_seq.replace('A','G')
        #print("after_editing_seq:\t", after_editing_seq)
    ### ref_full_seq, after_full_seq:
    ref_full_seq12=snp.seq5+snp.wt+snp.seq3## under_full_seq=snp.seq5+snp.mut+snp.seq3
    if under_editing_cor[0] != 0 and after_editing_seq:
        after_full_seq=snp.seq5[:under_editing_cor[0]]+after_editing_seq+snp.seq3[under_editing_cor[1]:] 
    elif under_editing_cor[0] == 0 and after_editing_seq:
        after_full_seq =snp.seq5+after_editing_seq+snp.seq3[under_editing_cor[1]:]

    if after_full_seq: 
        if rc_flag==False:
            aa_ref_full_seq=Seq(ref_full_seq12).translate()
            aa_after_full_seq=Seq(after_full_seq).translate()
        else: ### if rc_flag==True
            aa_ref_full_seq = Seq(ref_full_seq12).reverse_complement().translate()  
            aa_after_full_seq = Seq(after_full_seq).reverse_complement().translate()

        if aa_ref_full_seq==aa_after_full_seq:
            synonymous12_flag=True
        #print("ref_full_seq\t",ref_full_seq12,"\nlen(ref_full_seq)\t",len(ref_full_seq12))
        #print("after_full_seq\t", after_full_seq, "\nlen(after_full_seq)\t", len(after_full_seq))
        #print('synonymous12_flag:',synonymous12_flag)

    return synonymous12_flag
def precise_syn12_correct(snp):
    snp = cut_3nlen(snp) 
    #print(snp)
    ### ABE or CBE:
    rc_flag=False
    if snp.mut == "T":
        rc_flag=True
        BElist = ABElist
        rc_snp_seq5 = Seq(snp.seq3).reverse_complement()
        rc_snp_seq3 = Seq(snp.seq5).reverse_complement()
        rc_snp_mut="A"
        rc_snp_wt=Seq(snp.wt).reverse_complement()
        snp=snp_define(snp.jobid,rc_snp_seq5,rc_snp_seq3,rc_snp_wt,rc_snp_mut,snp.readingframe)
    elif snp.mut == "A":
        BElist = ABElist
    elif snp.mut == "G":
        rc_flag = True
        BElist = CBElist
        rc_snp_seq5 = Seq(snp.seq3).reverse_complement()
        rc_snp_seq3 = Seq(snp.seq5).reverse_complement()
        rc_snp_mut="C"
        rc_snp_wt=Seq(snp.wt).reverse_complement()
        snp=snp_define(snp.jobid,rc_snp_seq5,rc_snp_seq3,rc_snp_wt,rc_snp_mut,snp.readingframe)
    elif snp.mut == "C":
        BElist = CBElist
    else:
        BElist=None
    #######################################
    ### precise correction or synonymous12
    precise_snp_be_grna = {};synonymous12_snp_be_grna = {}
    precise_snp_be_grna[snp.jobid] = {};synonymous12_snp_be_grna[snp.jobid] = {}  
    precise_be_grna = {};synonymous12_be_grna = {}
    if BElist:
        for BE in BElist:  ## eg of BE: "BE4(max)": ["NGG",13,17,"C","T","U"]
            #print("BE:\n", BE)
            precise_be_grna[BE]=[];synonymous12_be_grna[BE] = []
            if BElist[BE][5] == "U":
                act_wind_range = BElist[BE][2] - BElist[BE][1] + 1
                pam_len = len(BElist[BE][0])
                BEpam = BElist[BE][0]
                ##################################
                #### if satisfy pam && basein_wind
                for i in range(-1, (act_wind_range - 2) + 1):
                    ## the most left:
                    # under_editing_seq=snp.seq5[-(act_wind_range-1):]+snp
                    # pam_seq=snp.seq3[BElist[BE][1]:(BElist[BE][1]+pam_len)]

                    #####  under_editing_seq, under_editing_cor, pam_seq, grna_seq
                    if (act_wind_range - i - 1 - 1) != 0:
                        under_editing_seq = snp.seq5[-(act_wind_range - i - 1 - 1):] + snp.mut + snp.seq3[:(i + 1)]
                    else:
                        under_editing_seq = snp.mut + snp.seq3[:(i + 1)]
                    under_editing_cor = [-(act_wind_range - i - 1 - 1), (i + 1)]  ## define snp_cor=0
                    pam_seq = snp.seq3[(BElist[BE][1] + i): (BElist[BE][1] + pam_len + i)]
                    if (20 - (BElist[BE][1] + i) - 1) != 0:
                        grna_seq = snp.seq5[-(20 - (BElist[BE][1] + i) - 1):] + snp.mut + snp.seq3[:(BElist[BE][1] + i)]
                    else:
                        grna_seq = snp.mut + snp.seq3[:(BElist[BE][1] + i)]
                    #print("under_editing_seq\tpam_seq\tgrna_seq\n", under_editing_seq, "\t", pam_seq, "\t", grna_seq)

                    ### if pam match:
                    pamFlag = False
                    base_count = 0
                    if len(pam_seq) == pam_len: 
                        for j in range(pam_len):
                            if pam_seq[j] in BEletter[BEpam[j]]:
                                base_count += 1
                        if base_count == pam_len:
                            pamFlag = True
                            # print(pam_seq)

                    ### precise: only 1 snp.mut 'C'/'A' IN wind && snp.wt=='T'/'G'
                    ### synonymous12: >1 snp.mut 'C'/'A' IN wind
                    preciseFlag = False; synonymous12Flag = False
                    if len(under_editing_seq.replace(BElist[BE][3], '')) + 1 == len(under_editing_seq) and snp.wt == \
                            BElist[BE][4]:
                        preciseFlag = True
                        # print(under_editing_seq)

                    if len(under_editing_seq.replace(BElist[BE][3], '')) + 1 < len(under_editing_seq) or snp.wt != \
                            BElist[BE][4] : 
                        ## synonymous case 1 ,case2
                        #print("try synonymous12\n")
                        synonymous12Flag = synonymous12(snp,rc_flag, under_editing_cor, BElist,BE)

                    ################# output dict ################
                    pam_seq=pam_seq.lower()
                    if pamFlag == True and preciseFlag == True:
                        precise_be_grna[BE].append(grna_seq + pam_seq)
                    if pamFlag == True and synonymous12Flag == True:
                        synonymous12_be_grna[BE].append(grna_seq + pam_seq)

        precise_snp_be_grna[snp.jobid].update(precise_be_grna)
        synonymous12_snp_be_grna[snp.jobid].update(synonymous12_be_grna)
        #print(precise_snp_be_grna)
        return precise_snp_be_grna, synonymous12_snp_be_grna
###############################################################
def suposnp_reffullseq_cut_3nlen(snp,supo_snp,ref_full_seq):
    f1 = lambda x, i: x%((x // 3 - 1) * 3 + (i-1))
    f2 = lambda x,i: (x // 3 - 1) * 3 + (3-i)

    if snp.readingframe==1:
        #if supo_snp.readingframe==1:
         #   ref_full_seq_3n=ref_full_seq[]
        if supo_snp.readingframe==2:
            supo_snp.seq5=supo_snp.seq5[f1(26,2):]
            supo_snp.seq3=supo_snp.seq3[:f2(24,2)]
            ref_full_seq=ref_full_seq[f1(26,2):(26+1+f2(24,2))]
        if supo_snp.readingframe == 3:
            supo_snp.seq5=supo_snp.seq5[f1(27,3):]
            supo_snp.seq3=supo_snp.seq3[:f2(23,3)]
            ref_full_seq=ref_full_seq[f1(27,3):(27+1+f2(23,3))]
    if snp.readingframe==2:
        if supo_snp.readingframe==1:
            supo_snp.seq5=supo_snp.seq5[f1(24,1):]
            supo_snp.seq3=supo_snp.seq3[:f2(26,1)]
            ref_full_seq=ref_full_seq[f1(24,1):(24+1+f2(26,1))]
        #if supo_snp.readingframe==2:
         #   ref_full_seq_3n=ref_full_seq[]
        if supo_snp.readingframe == 3:
            supo_snp.seq5=supo_snp.seq5[f1(26,3):]
            supo_snp.seq3=supo_snp.seq3[:f2(24,3)]
            ref_full_seq=ref_full_seq[f1(26,3):(26+1+f2(24,3))]
    if snp.readingframe == 3:
        if supo_snp.readingframe==1:
            supo_snp.seq5=supo_snp.seq5[f1(23,1):]
            supo_snp.seq3=supo_snp.seq3[:f2(27,1)]
            ref_full_seq=ref_full_seq[f1(23,1):(23+1+f2(27,1))]
        if supo_snp.readingframe==2:
            supo_snp.seq5 = supo_snp.seq5[f1(24, 2):]
            supo_snp.seq3 = supo_snp.seq3[:f2(26, 2)]
            ref_full_seq = ref_full_seq[f1(24, 2):(24+1+f2(26, 2))]
    return supo_snp,ref_full_seq

def match_pam_wind(snp,base2edit,BElist):
    BE_windCorGrnaPam_d={}
    for BE in BElist:
        BE_windCorGrnaPam_d[BE]=[]
        if BElist[BE][5] == "U":
            act_wind_range = BElist[BE][2] - BElist[BE][1] + 1
            pam_len = len(BElist[BE][0])
            BEpam = BElist[BE][0]

            for i in range(-1, (act_wind_range - 2) + 1):
                matchFlag=False;grna_seq=None;pam_seq=None
                ## the most left:
                # under_editing_seq=snp.seq5[-(act_wind_range-1):]+snp
                # pam_seq=snp.seq3[BElist[BE][1]:(BElist[BE][1]+pam_len)]
                if (act_wind_range - i - 1 - 1) != 0:
                    under_editing_seq = snp.seq5[-(act_wind_range - i - 1 - 1):] + base2edit + snp.seq3[:(i + 1)]
                else:
                    under_editing_seq = base2edit + snp.seq3[:(i + 1)]
                under_editing_cor = [-(act_wind_range - i - 1 - 1), (i + 1)]  ## define snp_cor=0
                pam_seq = snp.seq3[(BElist[BE][1] + i): (BElist[BE][1] + pam_len + i)]
                if (20 - (BElist[BE][1] + i) - 1) != 0:
                    grna_seq = snp.seq5[-(20 - (BElist[BE][1] + i) - 1):] + base2edit + snp.seq3[:(BElist[BE][1] + i)]
                else:
                    grna_seq = base2edit + snp.seq3[:(BElist[BE][1] + i)]
                # print("under_editing_seq\tpam_seq\tgrna_seq\n", under_editing_seq, "\t", pam_seq, "\t", grna_seq)
            #### if satisfy
                ## if pam match
                pamFlag = False
                base_count = 0
                if len(pam_seq) == pam_len: ###!!!!!!!!
                    #print("pam_seq:", pam_seq, "\npam_len:", pam_len)
                    for j in range(pam_len):
                        if pam_seq[j] in BEletter[BEpam[j]]:
                            base_count += 1
                    if base_count == pam_len:
                        pamFlag = True
                else:
                    continue
                ## if base2edit in wind
                baseFlag=False
                #print("type of under_editing_seq",type(under_editing_seq))
                if len(under_editing_seq.replace(base2edit,''))<len(under_editing_seq):  
                    baseFlag=True
                if (pamFlag==True and baseFlag==True):
                    matchFlag = True
                if matchFlag and grna_seq and pam_seq: ### !!!!!!!!!!!!!!!!!
                    BE_windCorGrnaPam_d[BE].append([under_editing_cor,grna_seq,pam_seq])
    return BE_windCorGrnaPam_d

def syn3_correct(snp):
    ### supo_snp: jobid,seq5,seq3,wt,readingframe
    ### snp:      jobid,seq5,seq3,wt,mut,readingframe
    ###################################################
    ######## define: supo_snp1/2, ref_full_seq1/2 #############
    ### ref_Full_seq: cut_3n, rc,-> a variety--consecutive segment--used for translate()
    if snp.readingframe == 3:
        supo_snp1_seq5, supo_snp2_seq5 = snp.seq5[:-2],snp.seq5[:-1]
        supo_snp1_seq3, supo_snp2_seq3 = (snp.seq5[-1]+snp.mut+snp.seq3),(snp.mut+snp.seq3)
        supo_snp1_wt, supo_snp2_wt = snp.seq5[-2], snp.seq5[-1]
        supo_snp1_readingframe, supo_snp2_readingframe = 1,2
        ref_full_seq1=supo_snp1_seq5+supo_snp1_wt+(snp.seq5[-1]+snp.wt+snp.seq3);ref_full_seq2=supo_snp2_seq5+supo_snp2_wt+(snp.wt+snp.seq3)
    elif snp.readingframe == 2:
        supo_snp1_seq5, supo_snp2_seq5 = snp.seq5[:-1],(snp.seq5+snp.mut)
        supo_snp1_seq3, supo_snp2_seq3 = (snp.mut+snp.seq3),snp.seq3[1:]
        supo_snp1_wt, supo_snp2_wt = snp.seq5[-1],snp.seq3[0]
        supo_snp1_readingframe, supo_snp2_readingframe = 1,3
        ref_full_seq1 = supo_snp1_seq5 + supo_snp1_wt +(snp.mut+snp.seq3);ref_full_seq2=(snp.seq5+snp.wt)+supo_snp2_wt+supo_snp2_seq3
    elif snp.readingframe == 1:
        supo_snp1_seq5, supo_snp2_seq5 = (snp.seq5+snp.mut),(snp.seq5+snp.mut+snp.seq3[0])
        supo_snp1_seq3, supo_snp2_seq3 = snp.seq3[1:],snp.seq3[2:]
        supo_snp1_wt, supo_snp2_wt = snp.seq3[0],snp.seq3[1]
        supo_snp1_readingframe, supo_snp2_readingframe = 2,3
        ref_full_seq1 =(snp.seq5+snp.wt)+supo_snp1_wt+supo_snp1_seq3;ref_full_seq2=(snp.seq5+snp.wt+snp.seq3[0])+supo_snp2_wt+supo_snp2_seq3

    supo_snp_l=None;ref_full_seq_l=None
    try:
        supo_snp1=supo_snp_define(snp.jobid,supo_snp1_seq5,supo_snp1_seq3,supo_snp1_wt,supo_snp1_readingframe)
        supo_snp2 = supo_snp_define(snp.jobid, supo_snp2_seq5, supo_snp2_seq3, supo_snp2_wt, supo_snp2_readingframe)
        supo_snp_l=[supo_snp1,supo_snp2]
        ref_full_seq_l=[ref_full_seq1,ref_full_seq2]
    except:
        print("a wrong snpin for synonymous3")
    #############################################
    ######  {} synonymous3_snp_be_grna   #########
    synonymous3_snp_be_grna = {}
    synonymous3_be_grna = {}
    synonymous3_snp_be_grna[snp.jobid] = {}  

    for i in [0,1]:
        supo_snp,ref_full_seq=suposnp_reffullseq_cut_3nlen(snp,supo_snp_l[i],ref_full_seq_l[i]) 
        #print("i of supo_snp,ref_full_seq:",i,"\n***",supo_snp,"\n***",ref_full_seq)
        ###### BElist=ABE or CBE; if rc_flag: rc-supo_snp, rc-ref_full_seq
        rc_flag=False 
        if supo_snp.wt=='C':
            BElist=CBElist
        elif supo_snp.wt=='G':
            BElist = CBElist
            rc_snp_seq5 = Seq(supo_snp.seq3).reverse_complement()
            rc_snp_seq3 = Seq(supo_snp.seq5).reverse_complement()
            rc_snp_wt = "C"
            supo_snp = supo_snp_define(supo_snp.jobid, rc_snp_seq5, rc_snp_seq3, rc_snp_wt, supo_snp.readingframe)
            ref_full_seq=Seq(ref_full_seq).reverse_complement() 
            rc_flag=True
        elif supo_snp.wt == 'A':
            BElist = ABElist
        elif supo_snp.wt == 'T':
            BElist = ABElist
            rc_snp_seq5 = Seq(supo_snp.seq3).reverse_complement()
            rc_snp_seq3 = Seq(supo_snp.seq5).reverse_complement()
            rc_snp_wt = "A"
            supo_snp = supo_snp_define(supo_snp.jobid, rc_snp_seq5, rc_snp_seq3, rc_snp_wt, supo_snp.readingframe)
            ref_full_seq = Seq(ref_full_seq).reverse_complement() 
            rc_flag = True
        else:
            BElist=None
        if i==0 and BElist!=None:
            for BE in BElist:
                synonymous3_be_grna[BE]=[]
        elif i==1 and BElist!=None:
            for BE in BElist:
                if BE not in synonymous3_be_grna.keys():
                    synonymous3_be_grna[BE]=[]
        ###########################################
        ##### edit_target: supo_snp.wt
        ### condition:  if aa_after_full_seq(supo_snp.edited, mutbase.maybeEdited)==aa_ref_full_seq(supo_snp.wt, mutbase.wt)
        #try:
        #print(type(supo_snp))
        #print("supo_snp.wt***\t" ,supo_snp.wt,"BElist***\t",BElist)
        ######################################################
        ## match_pam_wind(): in: snp, base2edit, BElist :
        BE_windCorGrnaPam_d=match_pam_wind(supo_snp, supo_snp.wt, BElist)
        #print("BE_windCorGrnaPam_d:\t",BE_windCorGrnaPam_d)

        for BE in BE_windCorGrnaPam_d.keys():
            #synonymous3_be_grna[BE]=[] ^
            #print("BE:\n", BE)
            for windCorGrnaPam_l in BE_windCorGrnaPam_d[BE]: 
                under_editing_cor= windCorGrnaPam_l[0] 
                synonymous3_flag=False
                ### under_editing_seq, after_editing_seq;  ref_full_seq, after_full_seq:
                if under_editing_cor[0] != 0:
                    under_editing_seq = supo_snp.seq5[under_editing_cor[0]:] + supo_snp.wt + supo_snp.seq3[:under_editing_cor[1]] # + snp.mut +
                else:
                    under_editing_seq = supo_snp.wt + supo_snp.seq3[:under_editing_cor[1]]

                after_editing_seq = None;after_full_seq = None  
                if BElist[BE][3] == 'C':
                    after_editing_seq = under_editing_seq.replace('C', 'T')
                    #print("after_editing_seq:\t", after_editing_seq)
                if BElist[BE][3] == 'A':
                    after_editing_seq = under_editing_seq.replace('A', 'G')
                    #print("after_editing_seq:\t", after_editing_seq)
                ### ref_full_seq, after_full_seq:
                if under_editing_cor[0] != 0 and after_editing_seq:
                    after_full_seq = supo_snp.seq5[:under_editing_cor[0]] + after_editing_seq + supo_snp.seq3[under_editing_cor[1]:]  
                elif under_editing_cor[0] == 0 and after_editing_seq:
                    after_full_seq = supo_snp.seq5 + after_editing_seq + supo_snp.seq3[under_editing_cor[1]:]

                if after_full_seq:  
                    #print("ref_full_seq:\t",ref_full_seq,"\nafter_full_seq:\t",after_full_seq)
                    if rc_flag==False:
                        aa_ref_full_seq = Seq(ref_full_seq).translate()  
                        aa_after_full_seq = Seq(after_full_seq).translate()
                    else:  ##ie rc_flag==True
                        aa_ref_full_seq = Seq(ref_full_seq).reverse_complement().translate()  
                        aa_after_full_seq = Seq(after_full_seq).reverse_complement().translate()
                    if aa_ref_full_seq == aa_after_full_seq:
                        synonymous3_flag = True
                        grna_seq=windCorGrnaPam_l[1]
                        pam_seq=windCorGrnaPam_l[2]
                        pam_seq = pam_seq.lower()
                        synonymous3_be_grna[BE].append(grna_seq + pam_seq)
            #print('synonymous3_be_grna:\t',synonymous3_be_grna) ## for i in [0,1]; for BE in BE_windCorGrnaPam_d
        synonymous3_snp_be_grna[snp.jobid].update(synonymous3_be_grna) ## for i in [0,1]
    #print('synonymous3_snp_be_grna:\t',synonymous3_snp_be_grna)
        #except:
         #   print("don't run synonymous3 search")
    return synonymous3_snp_be_grna
###################################################################
def dics_integrate(dic1,dic2):
    import copy
    dicout = copy.deepcopy(dic1)
    for snp in dicout:
        for be in dicout[snp]:
            if be in dic2[snp]:
                for grna in dic2[snp][be]:
                    dicout[snp][be].append(grna) if grna not in dicout[snp][be] else dicout[snp][be]
        for be in dic2[snp]:
            if be not in dicout[snp]:
                dicout[snp][be] = dic2[snp][be]
    return dicout
def re_extract(dic_snp_be):
    target_str=str(dic_snp_be)
    pattern=re.compile(r'\'([AGCTagct]{20,})\'')
    result=pattern.findall(target_str)
    return result
def filewrite(snpin,outPath,seq5,seq3, mut):
    precise = {};synonymous12 = {}
    synonymous3 = {}
    synonymous123={}

    synonymous3_snp_be_grna = syn3_correct(snpin)
    precise_snp_be_grna, synonymous12_snp_be_grna=precise_syn12_correct(snpin)
    print("********* precise_snp_be_grna *******",precise_snp_be_grna)
    precise[snpin.jobid]={}
    for be in precise_snp_be_grna[snpin.jobid]:
        if precise_snp_be_grna[snpin.jobid][be]!=[]:
            precise[snpin.jobid][be] =precise_snp_be_grna[snpin.jobid][be]

    synonymous12[snpin.jobid]={}
    for be in synonymous12_snp_be_grna[snpin.jobid]:
        if synonymous12_snp_be_grna[snpin.jobid][be]!=[]:
            synonymous12[snpin.jobid][be]=synonymous12_snp_be_grna[snpin.jobid][be]
    synonymous3[snpin.jobid] = {}
    for be in synonymous3_snp_be_grna[snpin.jobid]:
        if synonymous3_snp_be_grna[snpin.jobid][be] != []:
            synonymous3[snpin.jobid][be] = synonymous3_snp_be_grna[snpin.jobid][be]

    print("*****####### precise00000000000000:\n", precise)
    synonymous123=dics_integrate(synonymous12,synonymous3)
    allmatch=dics_integrate(precise,synonymous123)
    print("*****####### precise:\n", precise)
    print("*****####### synonymous12:\n", synonymous12)
    print("*****####### synonymous3:\n", synonymous3)
    print("********* synonymous123:\n",synonymous123)
    print("********* allmatch:\n", allmatch)
    ################  writeout  ##################
    with open("{0}/{1}_cleanMatches_file.txt".format(outPath,jobid), mode='w') as fout: 
        fout.write("snp\tprecise_BE-grna\n")
        #print("********precise*******",precise)
        dic=copy.deepcopy(precise) 
        #print("********dic of precise*******", dic)
        for snp in dic:
            fout.write("%s\t" % snp)
            for BE in dic[snp]:
                if dic[snp][BE] != []:
                    fout.write("%s:" % BE)
                    grna_l = re_extract(dic[snp][BE])
                    for grna in grna_l:
                        hive=hive_score(BE, seq5,seq3,mut, grna)
                    #print('grna_l:', grna_l)
                        fout.write("%s efficiency: %s " % (grna,hive))
                fout.write("\t")
            fout.write("\n")

    with open("{0}/{1}_quietMatches_file.txt".format(outPath,jobid), mode='w') as fout: 
        fout.write("snp\tsynonymous_BE-grna\n")
        dic=copy.deepcopy(synonymous123)
        for snp in dic:
            fout.write("%s\t" % snp)
            for BE in dic[snp]:
                if dic[snp][BE] != []:
                    fout.write("%s:" % BE)
                    grna_l = re_extract(dic[snp][BE])
                    for grna in grna_l:
                        hive = hive_score(BE, seq5,seq3,mut, grna)
                    #print('grna_l:', grna_l)
                        fout.write("%s efficiency: %s " % (grna, hive))
                fout.write("\t")
            fout.write("\n")


def BEsingle_100nt(jobid,seq5,seq3,wt,mut,readingframe,outPath):
    snpin = snp_define(jobid, seq5, seq3, wt, mut, readingframe)
    filewrite(snpin, outPath,seq5,seq3,mut)

def BEsingle_rsID(rsID,outPath):
    allsnps=get_rsIDsnp(rsID,outPath)
    try:
        for snp in allsnps:
            jobid='rs'+str(snp['formal_id'])
            seq5=snp['genomic_flanking'][0:50]
            seq3=snp['genomic_flanking'][51:]
            readingframe=int(snp['reading_frame'])
            if re.match(r'(\S{1})>(\S{1})',snp['SNP']):
                m=re.match(r'(\S{1})>(\S{1})',snp['SNP'])
                wt=m.group(1)
                mut=m.group(2)
            if wt and mut:
                snpin = snp_define(jobid, seq5, seq3, wt, mut, readingframe)
                filewrite(snpin, outPath,seq5,seq3,mut)

    except:
        print("Can't search SNV information of rsID:{0}".format(rsID))



"""
jobid, seq5, seq3, wt, mut, readingframe="rs80357410","GCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAG",\
                                         "GTCCCATCTGGTAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATG","T","C",1

outPath=""
snpin = snp_define(jobid, seq5, seq3, wt, mut, readingframe)
filewrite(snpin, outPath,seq5,seq3, mut)
"""
