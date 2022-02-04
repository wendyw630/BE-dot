#AUTHOR: Zelu Wang , wangzl630@163.com
from Bio.Seq import Seq
import csv
import sys,argparse
from Scripts.BEMain import BEMain
from Scripts.BEsingle import BEsingle
#    args
#    提供rsID,1.数据集有的，已经分析得到的clean quiet matches 结果文件
#            2.数据集没有的，提示输入50mer
#    ***提供50mer(upSeq, downSeq, mutation, wt, readingFrame), 参考脚本: besingle.py
#    允许用户custormize, BE类型(personalPAM,start,end,fromto,stream)，参考脚本：besingle.py, (BElist.update())
#    default=None

#    BEminorlist忽略

#    把rs 和50 mer两种输入形式分开写


def importCustomizeBE(CustomizeBE):
    with open(CustomizeBE,"r") as fin:
        for line in fin:
            line=line.strip()
            try:
                [BE_name,pam,start,end,BE_type,stream]=line.split(',')
                start=int(start);end=int(end)

                if BE_type=="ABE":
                    ABElist[BE_name]=[pam,start,end,"A","G", stream]
                elif BE_type=="CBE":
                    CBElist[BE_name] = [pam, start, end, "A", "G", stream]
            except:
                print("the input"+CustomizeBE+"is in a wrong format")
    return ABElist,CBElist

def if_snp_in_DB(DBfile,jobID,seq5,seq3,mut,wt,readingFrame):#in ".csv" format
    #requirment1. 51mer requirment2. snp-type

    ###########    read DB    ###########
    with open(DBfile,"r",encoding='UTF-8') as fin:
        snp51mer_dict={}
        csv_read=csv.reader(fin,delimiter=",")
        for line in csv_read:
            snp51mer_dict[line[7]]=line

    ###########   align 51mer with input seq5wtseq3, mut  ###########
    input_51mer=seq5+wt+seq3
    input_snp=wt+">"+mut
    if input_51mer in snp51mer_dict.keys() and input_snp==snp51mer_dict[input_51mer][9]:
        matchedlt=snp51mer_dict[input_51mer]
        fout=open(r'/home/wzl/test/predict/script/BE-FF/Scripts/snp_in_DB.csv',"w",newline='')
        csv_write = csv.writer(fout, dialect='excel') #dialect='excel'...........
        csv_write.writerow(matchedlt)
        fout.close()
        BEMain(r'/home/wzl/test/predict/script/BE-FF/Scripts/snp_in_DB.csv') #...........
    else:
        BEsingle(jobID,seq5,seq3,mut,wt,readingFrame) ## vary-version of BEMain

# ['\ufeff', 'SNP ID', 'Genome assembly', 'chromosome', 'orientation', 'start_pos', 'end_pos', '5Xmer', 'Clinical significance', 'SNP', 'residue', 'var type', 'Gene name', 'Gene ID', 'readingFrame', 'aaPosition', 'condition']
# ['rs28937596', 'rs28937596', 'GRCh38.p7', '3', '+', '184144111', '184144111', 'CTTTCTTCCATAGCTGCTAAAGGCCTGGAGCCCTGTTTTTAGGAACTACAT', 'pathogenic', 'T>C', 'W>R', 'missense', 'EIF2B5', '8893', '1', '628', 'Leukoencephalopathy with vanishing white matter (VWM)']


# get guide_RNA, PAM
def get_parser():
    parser = argparse.ArgumentParser(description="search available BE for a given snp, customized BE type is supported" )
    mer_input = parser.add_argument_group("# input snp context-50mer parameters")
    mer_input.add_argument("--jobID", metavar="<seq>", type=str,
                           help="give an ID to the search OTs job")
    mer_input.add_argument("--upSeq", metavar="<seq>", type=str,
                            help="25mer upstream of the snp")
    mer_input.add_argument("--downSeq", metavar="<seq>", type=str,
                           help="25mer downstream of the snp")
    mer_input.add_argument("--mutation", metavar="<seq>", type=str,
                           help="the mutation base")
    mer_input.add_argument("--wt", metavar="<seq>", type=str,
                           help="the ref base")
    mer_input.add_argument("--readingFrame", metavar="<int>", type=int,
                           help="the snp in codon readingFrame")


    customizeBE_input=parser.add_argument_group("# customized BE to add in BElist")
    customizeBE_input.add_argument("--CustomizeBE", metavar="<file>", type=str,default=None)
    #print(args.grna_input,args.pam_input,args.target_genome)
    return parser


def main():
    print("an eg. of snp 50mer input:")
    # upSeq, downSeq, mutation, wt, readingFrame,personalPAM,start,end,fromto,stream
    print("python3 BEsingle_50mer.py "
          " --jobID jobaaaaaaa"
          " --upSeq CCTGGACCCTAATGATTTTGACTTC"
          " --downSeq CGGTAACTGGGAGAGGGAGCCCCTC"
          " --mutation G --wt A --readingFrame 1")
    '''
    print("an eg. of snpID input:")
    print("python3 BEsingle_snpID.py "
          "--snpID rs61748411")
    '''


    print("an eg. of customized BE input:")
    print("python3 BEsingle.py --customizeBE personalBE.txt")
    ####### seeing example data:  ########
    # #BE_name,pam,start,end,BE_type,stream
    # BE4(max),NGG,13,17,CBE,U
    # BE3,NGNNR,13,18,CBE,U


    sys.stdout.write('#START: PIPELINE RUN HAS STARTED.\n')

    # STEP 1: Get the necessary arguments
    parser = get_parser()
    args = parser.parse_args()

    sys.stdout.write('#STEP 1: arguments parsed\n')
    for k in sorted(args.__dict__):
        if (args.__dict__[k] is not None) and args.__dict__[k] != '':
            sys.stdout.write('#ARG: ' + k + ' = ' + str(args.__dict__[k]) + '\n')

    # STEP 2:
    if args.CustomizeBE !=None:
        CustomizeBE=args.CustomizeBE
        importCustomizeBE(CustomizeBE)
    DBfile=r'/home/wzl/test/predict/script/BE-FF/snps_with_clinvar2.csv'
    [jobID, seq5, seq3, mut, wt, readingFrame]=[args.jobID,args.upSeq,args.downSeq,args.mutation,args.wt,args.readingFrame]
    if_snp_in_DB(DBfile, jobID, seq5, seq3, mut, wt, readingFrame)

if __name__=='__main__':
    main()


