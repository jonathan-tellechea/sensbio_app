import numpy as np
import pandas as pd
from rdkit import Chem

def seq_cod (seq):
    
    aa_list = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R",
                   "S","T","V","W","Y"]
    aam = np.zeros((20, 978)).astype(int)
    pos = 0
    for aa in seq:
        if aa == "*" or aa == "X":
            pos = pos + 1
            continue
        aa_idx = aa_list.index(aa)
        aam[aa_idx][pos] = 1
        pos = pos + 1
    aam = aam.reshape((1,20,978))
    
    return aam


def smi_cod(smi):
    
    ms = Chem.MolFromSmiles(smi)
    fp = Chem.RDKFingerprint(ms, fpSize = 512).ToBitString()
    fp = np.array([int(x) for x in list(fp)])
    fp = fp.reshape((1,512))
    
    return fp


def aam3d(db):
    
    aamatrix = seq_cod(db['AA_sequence'][0])
    count = 0
    for i in db['AA_sequence'][1:]:
        aam = seq_cod(i)
        aamatrix = np.concatenate((aamatrix, aam))
        count += 1
        print(count)
        
    return aamatrix


def fp3d(db):
    
    fpmatrix = smi_cod(db['SMILES'][0])
    count = 0
    for i in db['SMILES'][1:]:
        fpm = smi_cod(i)
        fpmatrix = np.concatenate((fpmatrix, fpm))
        count += 1
        print(count)
        
    return  fpmatrix

def pre_calc(flag):
    if flag == True:
        df_mol = pd.read_csv('./TF_DB_clean_pathway.csv')
        
        aamatrix = aam3d(df_mol)
        np.save('sequences_matrix.npy', aamatrix)
        
        fpmatrix = fp3d(df_mol)
        np.save('fingerprints_matrix.npy', fpmatrix)