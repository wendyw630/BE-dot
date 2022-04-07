#!/usr/bin/python
######### author: wangzelu ###########
######################################
from time import ctime
import argparse

from designsgRNA.utils import BEsingle_100nt
from designsgRNA.utils import BEsingle_rsID
from designsgRNA.baseEditorsTable import CBElist,ABElist,CGBElist
from OTprediction.utils import OT_predict
from OTannotation.utils import OT_anno
def init():
    desc = """Program: BE-dot
              Version: 1.0
              Author: Zelu Wang
              Email: wangzl630@163.com
    	"""
    usage = "%(prog)s"
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=desc,
                                     usage=usage)

    subparser = parser.add_subparsers(help="command", dest="command")

    design1parser = subparser.add_parser("designsgRNA_opt1",
                                        help="Design sgRNA for target SNV. \n  "
                                             "#Option 1: input the SNV and "
                                             "it's 100nt flanking sequence to design sgRNA")
    design1parser.add_argument('--jobID', required=True,type=str, help="Give an ID to the sgRNA-design job")
    design1parser.add_argument('--upSeq', required=True,metavar="<seq>", type=str, help="50nt upstream to the SNV")
    design1parser.add_argument('--downSeq', required=True,metavar="<seq>", type=str,help="50nt downstream to the SNV")
    design1parser.add_argument('--mut', required=True,type=str,help="Mutation base")
    design1parser.add_argument('--wt', required=True,type=str,help="Reference base")
    design1parser.add_argument('--codon_frame', metavar="<int>", type=int,
                              help="Amino acid codon frame that SNV occupy, 1/2/3")
    design1parser.add_argument('-o', '--outputPath',default="", type=str,
                               help="Path of output file(s)")


    design2parser = subparser.add_parser("designsgRNA_opt2",help="Design sgRNA for target SNV. \n  "
                                                                 "#Option 2: input the SNV's rsID to design sgRNA")
    design2parser.add_argument('--rsID', required=True,type=str, help="SNV's ID in dbSNP")
    design2parser.add_argument('-o', '--outputPath', default="", type=str,
                               help="Path of output file(s)")
    ###################  OT sites  ###########################
    BElist=[]
    for be in CBElist:
        BElist.append(be)
    for be in ABElist:
        BElist.append(be)
    for be in CGBElist:
        BElist.append(be)
    OTparser= subparser.add_parser("OTprediction",help="predict off-target profile")
    OTparser.add_argument('-BE',choices=BElist,required=True,help="BE to do OT-prediction")
    OTparser.add_argument('-grna',required=True,help="Single gRNA sequence to analyse (20nt)")
    #OTparser.add_argument('-onpam', required=True, help="pam of on-target site")
    OTparser.add_argument('-genome',"--target_genome", metavar="<file>", type=str,
                        help="Genome to search off-target sites",
                        default="Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa")
    OTparser.add_argument('-mis','--mismatch_number',default=3,
                          help="Maximum mismatch number considered in aligning. [default:%(default)s]")
    OTparser.add_argument('--DNAbulge',default=1,
                          help="Maximum DNA bulge number considered in aligning. [default:%(default)s]")
    OTparser.add_argument('--RNAbulge', default=1,
                          help="Maximum RNA bulge number considered in aligning. [default:%(default)s]")
    OTparser.add_argument('-o', '--outputPath', default="", type=str,
                               help="Path of output file(s)")
    ###################  ANNO  ###########################
    ANNOparser= subparser.add_parser("OTannotation",help="Annotate off-target products")
    ANNOparser.add_argument('-BE', "--base_editor", choices=BElist,help="BE to do OT-annotation")
    ANNOparser.add_argument('-i', "--input-file", help="Set input tsv file")
    ANNOparser.add_argument('-o', '--outputPath', default="", type=str,
                               help="Path of output file(s)")
    return parser

def main_cmd(args):
    now = datetime.now()
    if args.command == "designsgRNA_opt1":
        print("{0:40}{1:>20}".format("Design primer start!", ctime()))
        BEsingle_100nt(args.jobID,args.upSeq,args.downSeq,args.mut,args.wt,args.codon_frame,args.o) #
    if args.command == "designsgRNA_opt2":
        print("{0:40}{1:>20}".format("Design primer start!", ctime()))
        BEsingle_rsID(args.rsID,args.o) #
    if args.command == "OTprediction":
        print("{0:40}{1:>20}".format("Predict OT profile start!", ctime()))
        OT_predict(args.BE,args.grna,args.genome,args.mis,args.DNAbulge,args.RNAbulge,args.o)
    if args.command == "OTannotation":
        print("{0:40}{1:>20}".format("Annotate OT product start!", ctime()))
        OT_anno(args.BE,args.i,args.o)

if __name__ == "__main__":
    parser = init()
    args = parser.parse_args()

    if args.command not in ("designsgRNA_opt1","designsgRNA_opt2","formatdb","searchdb"):
        parser.print_help()
    else:
        main_cmd(args)
