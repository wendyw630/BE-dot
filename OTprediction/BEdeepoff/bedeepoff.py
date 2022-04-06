import argparse
import io
from pathlib import Path
import pickle

import re
import random
import numpy as np
import pandas as pd
import math
from collections import defaultdict, Counter


import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader,Sampler,BatchSampler,SubsetRandomSampler
import torch.nn.utils.rnn as rnn_utils

from Bio.Seq import Seq 
from Bio import pairwise2
from Bio.pairwise2 import format_alignment 
import utils

#Set computing environment
use_cuda = torch.cuda.is_available()
device = torch.device("cuda:0" if use_cuda else "cpu")

torch.set_printoptions(precision=6,sci_mode=False)
pd.set_option('display.float_format',lambda x : '%.6f' % x)

SEED = 1234

random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.cuda.manual_seed(SEED)

import argparse
parser = argparse.ArgumentParser(description="BEdeepoff for ABEmax and AncBE4max")
#parser.add_argument("-o","--off-output",help="Input Cas-OFFinder output or CRISPRitz output file",action='store_true')
parser.add_argument("-b", "--base-editor", 
                    choices=["ABE", "CBE"],
                    required=True,
                    help="Set base editor model")
parser.add_argument("-i","--input-file",help="Set input tsv file")
args = parser.parse_args()


stoi = {'<pad>':0, 'a':1,  'c':2,  'g':3,  't':4, 'A':1, 'C':2, 'G':3, 'T':4, '-':5}
itos = {0: '<pad>', 1:'a',  2:'c',  3:'g',  4:'t', 1:'A', 2:'C', 3:'G', 4:'T', 5:'-'}



print("Prepare computing environment...")
#Model configuration

ENC_INPUT_DIM = 6
ENC_EMB_DIM = 256
ENC_HID_DIM = 512
N_LAYERS = 1
ENC_DROPOUT = 0.5


class Net(nn.Module):
    def __init__(self, input_dim, emb_dim, enc_hid_dim, n_layers, dropout):
        super().__init__()
        
        self.enc_hid_dim = enc_hid_dim
        
        self.embedding = nn.EmbeddingBag(input_dim, emb_dim)
        
        self.rnn = nn.LSTM(emb_dim, enc_hid_dim, n_layers, bidirectional = True)
        
        
        self.dropout = nn.Dropout(dropout)
        self.fc_feat = nn.Linear(6*self.enc_hid_dim, 3*self.enc_hid_dim)
        self.fc_out = nn.Linear(3*self.enc_hid_dim, 1)
        self.len = 10
        
    def forward(self, seq, offset,length):
        global debug_var
        emb = self.dropout(self.embedding(seq, offset))
        
        #根据开始位置确定结束位置
        l_b = length.cumsum(dim=0)

        max_len = length[0]
        start_tmp = 0
        emb_tensors = []
        for i in range(len(length)):
            start = start_tmp
            end = l_b[i]
            emb_new = emb[start:end]
            seq_len = emb_new.shape[0]
            if max_len-seq_len>0:
                zeros = torch.zeros(max_len-seq_len,emb.shape[1]).to(device)
                emb_new = torch.cat([emb_new, zeros])
            emb_tensors.append(emb_new)
            start_tmp = end
            
        emb_tensor = torch.stack(emb_tensors).transpose(1,0)
        embed_input_x_packed = rnn_utils.pack_padded_sequence(emb_tensor, length, batch_first=False)

        pack_out, (hidden,_) = self.rnn(embed_input_x_packed)
        debug_var = pack_out,hidden, emb_tensor, embed_input_x_packed
        
        hidden = torch.cat( (hidden[-2,:,:], hidden[-1,:,:]), dim = 1 )
    
        pad_out, _ = rnn_utils.pad_packed_sequence(pack_out)
        avg_pool = torch.mean( pad_out.permute(1,0,2), 1)
        max_pool, _ = torch.max( pad_out.permute(1,0,2), 1)
        x = torch.cat([ hidden, avg_pool, max_pool ], dim=1)
        x = self.dropout(F.relu(self.fc_feat( x )))

        fc_out = self.fc_out(x)
        
        return fc_out

def get_model():
    abe_model_file = './Models/ABEdeepoff.pt'
    cbe_model_file = './Models/CBEdeepoff.pt'
    model1 = Net(ENC_INPUT_DIM, ENC_EMB_DIM, ENC_HID_DIM, N_LAYERS, ENC_DROPOUT).to(device)
    model1.load_state_dict( torch.load(abe_model_file,map_location=torch.device(device)))
    model2 = Net(ENC_INPUT_DIM, ENC_EMB_DIM, ENC_HID_DIM, N_LAYERS, ENC_DROPOUT).to(device)
    model2.load_state_dict( torch.load(cbe_model_file,map_location=torch.device(device)))
        
    return model1, model2


class gRNADataset(Dataset):
    def __init__(self,datafrm):
        self.df_data = datafrm
        df = self.df_data.apply(lambda x:get_encoding(x),axis=1)
        self.df = df.apply(pd.Series).reset_index(drop=True)
        self.df.columns = ['seq1','seq2','gRNA','target']
        self.indexes = self.df['seq1'].values
        self.offsets = self.df['seq2'].values
        self.gRNAs = self.df['gRNA']
        self.target = self.df['target']
    def __len__(self):
        
        return len(self.df)
    
    def __getitem__(self,idx):
        indexes = torch.LongTensor(self.indexes[idx])
        offsets = torch.LongTensor(self.offsets[idx])
        gRNA =  self.gRNAs[idx]
        target =  self.target[idx]
        return indexes, offsets, gRNA, target
    

#seq_idx, offset, y, a, c,off_type
def generate_batch(batch):

    indexes = []
    offsets = []
    seq_len = []
    old_max_ofs = 0
    gRNAs = []
    targets = []
    batch = [(a,b,b.shape[0],d, e)  for a, b,d,e in sorted(batch, key=lambda x: x[1].shape[0], reverse=True)]
    for i, (idx,ofs,lth,gRNA,target) in enumerate(batch):
        indexes += idx
        i = 1 if i > 0 else 0
        offset = old_max_ofs + ofs + i * 2
        old_max_ofs = offset[-1]
        offsets += list(offset)
        seq_len.append(lth)
        gRNAs.append(gRNA)
        targets.append(target)
        
    a = torch.LongTensor(indexes)
    b = torch.LongTensor(offsets)
    d = torch.LongTensor(seq_len)
    return a,b,d,gRNAs,targets  


def get_encoding(row):
    seq1 = row['gRNA']
    seq2 = row['target']
    alignments = pairwise2.align.globalms(seq1, seq2, 1, -1, -3, -2)
    a, b, c = format_alignment(*alignments[0]).split('\n')[:-2]
    seq_idx = []
    for i,(nuc1, alg, nuc2) in enumerate(zip(a,b,c)):
        seq_idx += [stoi[nuc1]] + [stoi[nuc2]]
    offset = [i for i in range(0,len(seq_idx),2)]
    return seq_idx, offset, a.replace('-',''), c.replace('-','')



def do_pred(iter_,model):
    model.eval()
    lst_dfs = []
    with torch.no_grad():
        for i, batch in enumerate(iter_):
            seq1 = batch[0].to(device)
            seq2 = batch[1].to(device)
            length = batch[2].to(device)

            out = model(seq1, seq2,length)
            out_eff = torch.sigmoid(out)
            out_eff = list(out_eff.view(-1).cpu().numpy())

            #获取gRNA 编辑结果的预测活性分布
            df_gRNA = pd.DataFrame({'source':batch[3],'target':batch[4]})
            df_gRNA['eff_pred'] = out_eff
            lst_dfs.append(df_gRNA)
    df_conc = pd.concat(lst_dfs)
    return df_conc

def prep_inputs(df_inputs):
    df_inputs['gRNA'] = df_inputs.source.str.upper().str.strip()
    df_inputs['target'] = df_inputs.target.str.upper().str.strip().str.replace('-','')
    df_inputs['target_len'] = df_inputs.target.apply(len)
    df_inputs.reset_index(drop=True,inplace=True)
    dataset = gRNADataset( df_inputs )
    batches =  DataLoader( dataset, batch_size=5, shuffle=False,
                      collate_fn=generate_batch )
    return batches

print("Load model...")
model_off = get_model()

def main():
    input_file = Path(args.input_file)
    stem = input_file.stem
    output_dir = utils.safe_makedir(input_file.parent / 'off_output')
    eff_file = input_file.parent / 'off_output' / (str(input_file.stem) + '_eff.csv')

    if not input_file.is_file():
        exit("File doesn't exist!")
    else:
        print("Load input file...")
        try:
            #CBE
            if input_file:
                df_inputs = pd.read_csv( input_file, sep='\t')
                if len( df_inputs.columns ) == 8:
                    '''
                    Cas-Offinder output
                    '''
                    cols = ['#Bulge_type', 'crRNA', 'DNA', 'Chromosome', 'Position', 'Direction',
        'Mismatches', 'Bulge_Size']
                    df_inputs.columns = cols                
                else:
                    '''
                    CRISPRtiz output
                    '''
                    cols = ['#Bulge_type', 'crRNA', 'DNA', 'Chromosome', 'Position',
        'Cluster Position', 'Direction', 'Mismatches', 'Bulge_Size', 'Total']
                    df_inputs.columns = cols
                    df_inputs = df_inputs[['#Bulge_type', 'crRNA', 'DNA', 'Chromosome', 'Position',
        'Direction', 'Mismatches', 'Bulge_Size']]
                    
                source = df_inputs[df_inputs.Mismatches == 0].DNA.values[0]
                df_inputs['source'] = source
                df_inputs['target'] = df_inputs.DNA
                df_inputs.reset_index(drop=True,inplace=True) 
            else:
                df_inputs = pd.read_csv( input_file, sep=' ',header=None )
                df_inputs.columns = ['source','target']
            batches = prep_inputs(df_inputs)
            model = model_off[0] if args.base_editor == "ABE" else model_off[1]
            print("Do prediction...")
            df_eff = do_pred( batches, model ).reset_index( drop=True )
            if input_file:
                df_eff = pd.concat([ df_inputs.iloc[:,:8], df_eff.eff_pred ], axis=1)
            df_eff.sort_values( by='eff_pred', ascending=False ).to_csv(eff_file)
            print("Finished prediciton!")
        except Exception as e:
            print(str(e))

if __name__ == "__main__":
    main()
