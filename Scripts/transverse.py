##this is the special case of checking the transverse bases
from Scripts.BEseqType import *
from Bio.Seq import Seq
from Scripts.baseEditorsTable import CBElist, CBElistMinor, ABElist, ABElistMinor, BEletter
#from Scripts.BEMain import cutFromRF,getRevComp

def printPamSeq(seq5,activationWindow,seq3, locfromEnd,PAM):
    len3=len(seq3)
    printPAM=seq5+"<class style='color:blue'><b>"+activationWindow+"</b></class>"+seq3[0:len3-locfromEnd]+"<b><class style='color:#DD96F0'>"+seq3[len3-locfromEnd:len3-locfromEnd+len(PAM)]+"</b></class>"+seq3[len3-locfromEnd+len(PAM):]
    return printPAM

def SpecialCleanMatch(snp, Matches, BElist,rev,originalProtein):
    cleanMatches={}
    clean_list=[]
    quietDic={}
    quiet_list=[]
    locations_dic={}
    origProtein_dic={}
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
            diff = 0


            PAM=BElist[BE][0]
            start=BElist[BE][1]-1
            end=BElist[BE][2]-1
            locations=[]  #find locations of PAM. for each location, check if clean
            window = (end - start+1)
            for i in range(window):
                temp = start+i
                lenPAM= len(PAM)
                j=0
                while (j<=lenPAM):
                    if seq3[temp] in BEletter[PAM[j]]:
                        j+=1
                        temp+=1
                    else:
                        break
                    if j == len(PAM):
                        locations.append(i + start) #maybe j is not neccesary?
                        break
            locations_dic[BE] = locations

            totalSeq1 = seq5 + snp.mutation + seq3
            totalSeq2 = seq5 + snp.wt + seq3
            printSeq5 = seq5  # for results
            printSeq3 = seq3  # for results
            diff = 0
            # if snp.readingFrame == 1:
            #     print(snp.snpID, "reading frame is 1")
            #     if len(snp.seq5) % 3 == 1:
            #         print ("take 1")
            #         totalSeq1 = totalSeq1[1:]
            #         totalSeq2 = totalSeq2[1:]
            #         printSeq5 = printSeq5[1:]
            #         print ("totalseq", totalSeq1)
            #         diff = 1
            #     elif len(snp.seq5) % 3 == 2:
            #         totalSeq1 = totalSeq1[2:]
            #         totalSeq2 = totalSeq2[2:]
            #         printSeq5 = printSeq5[2:]
            #         diff = 2
            # elif snp.readingFrame == 2:
            #     if len(snp.seq5) % 3 == 0:
            #         totalSeq1 = totalSeq1[2:]
            #         totalSeq2 = totalSeq2[2:]
            #         printSeq5 = printSeq5[2:]
            #         diff = 2
            #     elif len(snp.seq5) % 3 == 2:
            #         totalSeq1 = totalSeq1[1:]
            #         totalSeq2 = totalSeq1[1:]
            #         printSeq5 = printSeq5[1:]
            #         diff = 1
            # elif snp.readingFrame == 3:
            #     if len(snp.seq5) % 3 == 0:
            #         totalSeq1 = totalSeq1[1:]
            #         totalSeq2 = totalSeq2[1:]
            #         printSeq5 = printSeq5[1:]
            #         diff = 1
            #     elif len(snp.seq5) % 3 == 1:
            #         totalSeq1 = totalSeq1[2:]
            #         totalSeq2 = totalSeq2[2:]
            #         printSeq5 = printSeq5[2:]
            #         diff = 2
            # if len(totalSeq1) % 3 == 1:
            #     totalSeq1 = totalSeq1[:-1]
            #     totalSeq2 = totalSeq2[:-1]
            #     printSeq3 = printSeq3[:-1]
            # if len(totalSeq1) % 3 == 2:
            #     totalSeq1 = totalSeq1[:-2]
            #     totalSeq2 = totalSeq2[:-2]
            #     printSeq3 = printSeq3[:-2]
            if rev==False:
                protein_seq = Seq(totalSeq2).translate()
            else:
                protein_seq = Seq(totalSeq2).reverse_complement().translate()
            origProtein_dic[BE]=protein_seq
            #clean_list = []
            #quiet_list = []
            max_num = 0
            protein_match = False
            for loc in locations_dic[BE]:
                protein_match = False
                temp_seq=totalSeq1
                # loc=loc+len(snp.seq5)
                # activation_window = totalSeq1[loc-end-diff:loc-start]
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

                # finalSeq=totalSeq1[0:loc-end-diff]+new_AW+totalSeq1[loc-start:]
                # protein_seq_new=Seq(finalSeq).translate()

                if rev==False:
                    finalSeq = Seq(totalSeq1[0:loc - end] + str(new_AW) + totalSeq1[loc - start + 1:])
                else:
                    beginningP = Seq(
                        totalSeq1[0:len(printSeq3) - locFromEnd + len(printSeq5) - end]).reverse_complement()
                    new_AW = Seq(new_AW).reverse_complement()
                    # endP=Seq(totalSeq1[loc - start:]).reverse_complement()
                    endP = Seq(
                        totalSeq1[len(printSeq3) - locFromEnd + len(printSeq5) - start + 1:]).reverse_complement()
                    finalSeq = endP + new_AW + beginningP
                # if snp.readingFrame=="1":
                #     protein_seq_new = finalSeq[1:].translate()
                # if snp.readingFrame=="2":
                #     protein_seq_new = finalSeq.translate()
                # if snp.readingFrame=="3":
                protein_seq_new = finalSeq.translate()
                if max_num==1:
                    clean_list.append(BE)
                if originalProtein==protein_seq_new:
                    protein_match=True

                if protein_match==True:
                    quiet_list.append(BE)

    return clean_list,quiet_list


def origPro(snp,rev):
    totalSeq2 = snp.seq5 + snp.wt + snp.seq3
    if snp.readingFrame == "1":
        if len(snp.seq5) % 3 == 1:
            totalSeq2 = totalSeq2[1:]
        elif len(snp.seq5) % 3 == 2:
            totalSeq2 = totalSeq2[2:]
    elif snp.readingFrame == "2":
        if len(snp.seq5) % 3 == 0:
            totalSeq2 = totalSeq2[2:]
        elif len(snp.seq5) % 3 == 2:
            totalSeq2 = totalSeq2[2:]
    elif snp.readingFrame == "3":
        if len(snp.seq5) % 3 == 0:
            totalSeq2 = totalSeq2[1:]
        elif len(snp.seq5) % 3 == 1:
            totalSeq2 = totalSeq2[2:]

    if len(totalSeq2) % 3 == 1:
        totalSeq2 = totalSeq2[:-1]
    if len(totalSeq2) % 3 == 2:
        totalSeq2 = totalSeq2[:-2]
    totalSeq2 = Seq(totalSeq2)
    if rev == False:
        protein_seq = totalSeq2
    else:
        protein_seq = totalSeq2.reverse_complement()
    return protein_seq.translate()