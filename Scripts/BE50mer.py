#!/usr/bin/python3
from Scripts.BEseqType_single import *
from Bio.Seq import Seq
import csv
from Scripts.baseEditorsTable import CBElist, CBElistMinor, ABElist, ABElistMinor, BEletter
from Scripts.transverse import SpecialCleanMatch,origPro

#def importSNPS(SNPfile):

def import50merSNP(jobID,seq5,seq3,mut,wt,readingFrame):
    # rslts=[]
    newSNP=BEseqType(jobID,seq5,seq3,mut,wt,readingFrame)
    return newSNP

def matchBE(snp, BElist): #returns all the BE who's PAM'S are found in the given sequence
    matches_list = []
    matches = {}
    for BE in BElist:
        if BElist[BE][5] == "U":### wzl_note: if BElist[BE][5]=='U',sequence=snp.seq3; else: seq5[::-1]
            #################### wzl_note: only BElist[Cas12a-BE][5] == "D"
            sequence=snp.seq3
            #print("'BE' is:",BE, "'sequence' is snp.seq3:", sequence)  ###wzl_added
        else:
            sequence=snp.seq5
            sequence=sequence[::-1] # get reverse
        #print("'sequence' is ",sequence)###wzl_added

        PAM = BElist[BE][0]
        # find if PAM matches in correct place
        if len(sequence) > BElist[BE][2]:
            match = False
            start = BElist[BE][1] -1
            end = BElist[BE][2]-1
            window = (end-start+1)
            for i in range(window):
                temp = start+i
                lenPAM= len(PAM)
                j=0
                while (j<lenPAM):### wzl_changed, original: j<=lenPAM,IndexError: string index out of range
                    if sequence[temp] in BEletter[PAM[j]]:
                        j+=1
                        temp+=1
                    else:
                        break
                    if j == len(PAM): ### wzl_note: the condition fullfilled only if 3 consecutive-base string of sequence[12,13,14,15,16] is 'NGG'
                        match = True
                        break
            if match is True and BE!= "eA3A-BE3": ###wzl_note: a CBE, edit when C comes after T
                matches_list.append(BE)
            elif match is True and snp.seq5[-1]=='T': ### wzl_note: ie:
                                                      # (match is True) and (BE== "eA3A-BE3") and (snp.seq5[-1]=='T')
                matches_list.append(BE)
        matches[snp.jobID] = matches_list
    return matches

## This function takes as input the found matches_list, and checks for precise corrections
def cleanMatch(snp, Matches, BElist,rev):
    cleanMatches={}
    clean_list=[]
    quietDic={}
    quiet_list=[]
    locations_dic={}
    snp_cleanbe = {}  ### wzl_add
    snp_cleanbe[snp.jobID] = {}  ### wzl_add
    be_grna = {}  ### wzl_add
    for match in Matches:
        for BE in Matches[match]:

            if BElist[BE][5] == "U":
                seq3 = snp.seq3
                seq5=snp.seq5
            else:
                seq3 = snp.seq5[::-1]  # get reverse
                seq5 = snp.seq3[::-1]  # get reverse

            totalSeq1 = seq5 + snp.mutation + seq3
            totalSeq2 = seq5 + snp.wt + seq3
            printSeq5 = seq5  # for results
            printSeq3 = seq3  # for results

            if rev == False:
                protein_seq = Seq(totalSeq2).translate()
            else:
                protein_seq = Seq(totalSeq2).reverse_complement().translate()
            origProtein = protein_seq

            PAM=BElist[BE][0]
            start=BElist[BE][1]-1
            end=BElist[BE][2]-1
            locations=[]
            #find locations of PAM. for each location, check if clean
            window = (end - start+1)
            lenPAM = len(PAM)
            for i in range(window):
                temp = start+i

                j=0
                while (j<=lenPAM):
                    if seq3[temp] in BEletter[PAM[j]]:
                        j+=1
                        temp+=1
                    else:
                        break
                    if j == len(PAM):
                        locations.append(i + start)
                        break
            locations_dic[BE] = locations

            #offtarget={} ### wzl_add !!!!!!!!!!!!!!!!!!!!!
            offtarget_pam={} ### wzl_add !!!!!!!!!!!!!!!!!!!!!
            for loc in locations_dic[BE]:
                protein_match = False
                locFromEnd = len(printSeq3) - loc
                loc = loc + len(printSeq5)
                activation_window = totalSeq1[len(printSeq3) - locFromEnd + len(printSeq5) - end:len(
                    printSeq3) - locFromEnd + len(printSeq5) - start + 1]
                num = 0  # number of times the variant appears within the activation window
                max_num = 0
                new_AW = []
                for i in range(len(activation_window)):
                    if activation_window[i]==snp.mutation:
                        num=num+1
                        new_AW.append(snp.wt)
                    else:
                        new_AW.append(activation_window[i])
                new_AW=''.join(new_AW)
                if num>max_num:
                    max_num=num

                a = printSeq5 + snp.mutation + printSeq3 ### wzl_add
                pam_seq = a[loc + 1:loc + 1 + lenPAM].lower()  ### wzl_add, type: <str>
                offtarget_pam[loc] = a[loc - 19:loc + 1]+pam_seq ### wzl_add !!!!!!!!!!!!!!!!!!!!!


                if rev==False:
                    finalSeq = Seq(totalSeq1[0:loc - end] + str(new_AW) + totalSeq1[loc - start + 1:])
                else:
                    beginningP = Seq(
                        totalSeq1[0:len(printSeq3) - locFromEnd + len(printSeq5) - end]).reverse_complement()
                    new_AW = Seq(new_AW).reverse_complement()
                    endP = Seq(
                        totalSeq1[len(printSeq3) - locFromEnd + len(printSeq5) - start + 1:]).reverse_complement()
                    finalSeq = endP + new_AW + beginningP
                protein_seq_new = finalSeq.translate()
                if max_num==1:
                    clean_list.append(BE)
                    ##### wzl_note: clean_list  !!!!!!
                if protein_seq==protein_seq_new:
                    protein_match=True

                if protein_match==True:
                    quiet_list.append(BE)
                    ##### wzl_note: quiet_list  !!!!!!

            offtarget_pam_l=[] ### wzl_add
            for loc in offtarget_pam: ### wzl_add
                offtarget_pam_l.append(offtarget_pam[loc])
            be_grna[BE]=offtarget_pam_l ### wzl_add

        cleanbe_grna={}
        for clean_be in clean_list:
            cleanbe_grna[clean_be]=be_grna[clean_be]
        snp_cleanbe[snp.jobID].update(cleanbe_grna) ### wzl_add

        cleanMatches[snp.jobID]=snp_cleanbe[snp.jobID]### wzl_add  orig: cleanMatches[snp.jobID]=clean_list
        quietDic[snp.jobID]=quiet_list
    return cleanMatches, quietDic,origProtein

def beginningCut(snp):
    if snp.readingFrame == "1":
        if len(snp.seq5) % 3 == 1:
            snp.seq5 = snp.seq5[1:]
        elif len(snp.seq5) % 3 == 2:
            snp.seq5 = snp.seq5[2:]
    elif snp.readingFrame == "2":
        if len(snp.seq5) % 3 == 0:
            snp.seq5 = snp.seq5[2:]
        elif len(snp.seq5) % 3 == 2:
            snp.seq5 = snp.seq5[1:]
    elif snp.readingFrame == "3":
        if len(snp.seq5) % 3 == 0:
            snp.seq5 = snp.seq5[1:]
        elif len(snp.seq5) % 3 == 1:
            snp.seq5 = snp.seq5[2:]

    if (len(snp.seq5)+len(snp.seq3)+1) % 3 == 1:
        snp.seq3 = snp.seq3[:-1]
    if (len(snp.seq5)+len(snp.seq3)+1) % 3 == 2:
        snp.seq3 = snp.seq3[:-2]

    return snp

#def cutFromRF(snp,totalSeq1,totalSeq2,printSeq5,printSeq3)
### wzl_note: delete this function <- in original script, not been used
def getRevComp(snp):
    len3=len(snp.seq3)
    totalSeq= snp.seq5+snp.mutation+snp.seq3
    total_seq=Seq(totalSeq)
    rc_seq=total_seq.reverse_complement()
    seq5=str(rc_seq[0:len3])
    seq3=str(rc_seq[len3+1:])
    new_snp=BEseqType(snp.jobID,seq5 ,seq3,snp.mutation,snp.wt, snp.readingFrame)
    return new_snp

def checkRF(snp):
    #this function will check the other 2 bases in the reading frame to see whether fixing them may result in a quiet result
    if snp.readingFrame=="1":
        zero_seq5 = snp.seq5
        zero_mutation = snp.mutation
        zero_wt = find_cor(zero_mutation)
        zero_seq3 = snp.seq3
        snp0 = BEseqType(snp.jobID, zero_seq5, zero_seq3, zero_wt, zero_mutation, 1)
        first_seq5=snp.seq5+snp.mutation
        first_mutation=snp.seq3[0]
        first_wt=find_cor(first_mutation)
        first_seq3=snp.seq3[1:]
        snp1=BEseqType(snp.jobID,first_seq5,first_seq3,first_wt,first_mutation,2)
        second_seq5 = snp.seq5 + snp.mutation+snp.seq3[0]
        second_mutation = snp.seq3[1]
        second_wt = find_cor(second_mutation)
        second_seq3 = snp.seq3[2:]
        snp2 = BEseqType(snp.jobID, second_seq5, second_seq3, second_wt,second_mutation, 3)
        return snp0, snp1,snp2

    elif snp.readingFrame=="2":
        zero_seq5 = snp.seq5
        zero_mutation = snp.mutation
        zero_wt = find_cor(zero_mutation)
        zero_seq3 = snp.seq3
        snp0 = BEseqType(snp.jobID, zero_seq5, zero_seq3, zero_wt, zero_mutation, 2)
        first_seq5=snp.seq5[:-1]
        first_mutation=snp.seq5[-1]
        first_wt=find_cor(first_mutation)
        first_seq3=snp.mutation+snp.seq3
        snp1=BEseqType(snp.jobID,first_seq5,first_seq3,first_wt,first_mutation,1)
        second_seq5 = snp.seq5 + snp.mutation
        second_mutation = snp.seq3[0]
        second_wt = find_cor(second_mutation)
        second_seq3 = snp.seq3[1:]
        snp2 = BEseqType(snp.jobID, second_seq5, second_seq3, second_wt,second_mutation, 3)
        return snp0,snp1,snp2

    # elif snp.readingFrame==3:
    else:
        zero_seq5 = snp.seq5
        zero_mutation = snp.mutation
        zero_wt = find_cor(zero_mutation)
        zero_seq3 = snp.seq3
        snp0 = BEseqType(snp.jobID, zero_seq5, zero_seq3, zero_wt, zero_mutation, 3)
        first_seq5=snp.seq5[:-2]
        first_mutation=snp.seq5[-2]
        first_wt=find_cor(first_mutation)
        first_seq3=snp.seq5[-1]+snp.mutation+snp.seq3
        snp1=BEseqType(snp.jobID,first_seq5,first_seq3,first_wt,first_mutation,1)
        second_seq5 = snp.seq5[:-1]
        second_mutation = snp.seq5[-1]
        second_wt = find_cor(second_mutation)
        second_seq3 = snp.mutation+snp.seq3
        snp2 = BEseqType(snp.jobID, second_seq5, second_seq3,second_wt, second_mutation,2)
        return snp0, snp1,snp2

def find_cor(base):
    if base=="C":
        return "T"
    if base=="T":
        return "C"
    if base=="A":
        return "G"
    if base=="G":
        return "A"

def Main(jobID,seq5,seq3,mut,wt,readingFrame):
    snp=import50merSNP(jobID,seq5,seq3,mut,wt,readingFrame) #parsing cvs file
    matches={}
    minor_clean={}
    minor_quiet={}
    cleanMatchdic = {}
    quietMatchdic={}
    # sort DNA. First, determine which bases we wish to replace. 4 cases:
    # 1. C to T: use CBE list       2. A to G: use ABE list
    # 3. T to C: switch to reverse complement and use ABE
    # 4. G to A: switch to RC and use CBE

    #print("this is snp--",snp)### wzl_added
    rev = False
    snp = beginningCut(snp) ### wzl_note: snp.seq5, snp.seq3 are cut based on snp.readingFrame
    if snp.mutation == "C" and snp.wt == "T":
        BElist = CBElist
        MinorBElist=CBElistMinor
    elif snp.mutation == "A" and snp.wt == "G":
        BElist = ABElist
        MinorBElist = ABElistMinor
    elif snp.mutation == "T" and snp.wt == "C":
        snp=getRevComp(snp) ##### wzl_note: the contents of snp: seeing getRevComp(snp)
        rev = True
        snp.mutation = "A"
        snp.wt="G" ###wzl_note: 4 above steps, make snp from +strand's C->T to -strand's G->A
        BElist = ABElist
        MinorBElist = ABElistMinor
    elif snp.mutation == "G" and snp.wt == "A":
        snp = getRevComp(snp)
        rev = True
        snp.mutation = "C"
        snp.wt="T"
        BElist = CBElist
        MinorBElist = CBElistMinor
    else:
        BElist=None
        MinorBElist=None ###wzl_!!!

    try:
      # check for matches in major window
        check_match = matchBE(snp, BElist)
        #print('###### end of match BElist:') ###wzl_added
        matches.update(check_match)
    except:
        pass
    try:
      # check for matches in minor window
        check_match_minor = matchBE(snp, MinorBElist)
        #print('###### end of match_minor list')  ### wzl_added
        matches.update(check_match_minor)
    except:
        pass
    try:
        snp0, snp2, snp3 = checkRF(snp)
        rev0,rev2,rev3=False,False,False
        if snp0.mutation == "C":
            BElist0 = CBElist
        elif snp0.mutation == "A":
            BElist0 = ABElist
        elif snp0.mutation == "T":
            snp0 = getRevComp(snp0)
            rev0 = True
            snp0.mutation = "A"
            snp0.wt = "G"
            BElist0 = ABElist
        elif snp0.mutation == "G":
            snp0 = getRevComp(snp0)
            rev0 = True
            snp0.mutation = "C"
            snp0.wt = "T"
            BElist0 = CBElist
        if snp2.mutation == "C":
            BElist2 = CBElist
        elif snp2.mutation == "A":
            BElist2 = ABElist
        elif snp2.mutation == "T":
            snp2 = getRevComp(snp2)
            rev2 = True
            snp2.mutation = "A"
            snp2.wt = "G"
            BElist2 = ABElist
        elif snp2.mutation == "G":
            snp2 = getRevComp(snp2)
            rev2 = True
            snp2.mutation = "C"
            snp2.wt = "T"
            BElist2 = CBElist

        if snp3.mutation == "C":
            BElist3 = CBElist
        elif snp3.mutation == "A":
            BElist3 = ABElist
        elif snp3.mutation == "T":
            snp3 = getRevComp(snp3)
            rev3 = True
            snp3.mutation = "A"
            snp3.wt = "G"
            BElist3 = ABElist
        elif snp3.mutation == "G":
            snp3 = getRevComp(snp3)
            rev3 = True
            snp3.mutation = "C"
            snp3.wt = "T"
            BElist3 = CBElist
    except:
        pass
    try:
        check_match0 = matchBE(snp0, BElist0)
        matches.update(check_match0)
    except:
        pass
    try:
        check_match2 = matchBE(snp2, BElist2)
        matches.update(check_match2)
    except:
        pass
    try:
        check_match3 = matchBE(snp3, BElist3)
        matches.update(check_match3)
    except:
        pass

        #check for matches in minor window
        #check_match_minor = matchBE(snp, MinorBElist)
        #matchesMinor.update(check_match_minor)

    #check for clean match
    try:
        clean, quiet,orig_protein= cleanMatch(snp,check_match,BElist,rev)
    except:
        clean={}
        quiet={}
        orig_protein=origPro(snp,rev)

    try:
        clean_m, quiet_m, orig_protein_m = cleanMatch(snp, check_match_minor, MinorBElist, rev)
        for key in clean_m:
            minor_clean[key] = []
            for BE in clean_m[key]:
                if BE in clean[key]:
                    minor_clean[key].append(BE)
        for key in quiet_m:
            minor_quiet[key] = []
            for BE in quiet_m[key]:
                if BE in quiet[key]:
                    minor_quiet[key].append(BE)
    except:
        pass
    try:
        clean0,quiet0=SpecialCleanMatch(snp0,check_match0,BElist0,rev0,orig_protein)
    except:
        quiet0=[]
    try:
        clean2, quiet2 = SpecialCleanMatch(snp2, check_match2, BElist2,rev2,orig_protein)
    except:
        quiet2 = []
    try:
        clean3, quiet3 = SpecialCleanMatch(snp3, check_match3, BElist3,rev3,orig_protein)
    except:
        quiet3 = []

    for key in quiet0:
        if snp.jobID not in quiet:
            quiet[snp.jobID]=[key]
        elif key not in quiet[snp.jobID]:
            quiet[snp.jobID].append(key)
    for key in quiet2:
        if snp.jobID not in quiet:
            quiet[snp.jobID]=[key]
        elif key not in quiet[snp.jobID]:
            quiet[snp.jobID].append(key)
    for key in quiet3:
        if snp.jobID not in quiet:
            quiet[snp.jobID]=[key]
        elif key not in quiet[snp.jobID]:
            quiet[snp.jobID].append(key)
    try:
        quietMatchdic.update(quiet)
        cleanMatchdic.update(clean)
    except:
        pass
## statistics ##
## wzl_note: delete this part
    return matches, cleanMatchdic, quietMatchdic, minor_clean,minor_quiet
def BEsingle(jobID,seq5,seq3,mut,wt,readingFrame):
    matches, cleanMatchdic, quietMatchdic ,minor_clean,minor_quiet=Main(jobID,seq5,seq3,mut,wt,readingFrame)
    with open('single_matches_file.csv', mode='w') as f:
        f.write('jobID, matches\n')
        for key in matches.keys():
            f.write("%s,%s\n" % (key, matches[key]))

    with open('single_cleanMatches_file.csv', mode='w') as f:
        f.write('jobID, matches\n')
        #f.write(str(len(cleanMatchdic)))
        for rsid in cleanMatchdic.keys():
            f.write("%s\n"%rsid)
            for be in cleanMatchdic[rsid].keys():
                f.write("%s,%s\n" % (be, cleanMatchdic[rsid][be]))

    with open('single_quietMatches_file.csv', mode='w') as f:
        f.write('jobID, matches\n')
        #f.write(str(len(quietMatchdic)))
        for key in quietMatchdic.keys():
            f.write("%s,%s\n" % (key, quietMatchdic[key]))
