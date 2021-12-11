import streamlit as st
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def app():
    st.markdown("""
    # Molecular similarity tool

    Use this page to find similar entries in the database using your target compound in SMILES format as input.

    """)
    #defines a function to calculate the tanimoto similarity between two molecules in SMILES format (from https://medium.com/data-professor/how-to-calculate-molecular-similarity-25d543ea7f40).
    def tanimoto_calc(smi1, smi2):
        mol1 = Chem.MolFromSmiles(smi1)
        mol2 = Chem.MolFromSmiles(smi2)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 5, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 5, nBits=2048)
        sim = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
        return sim

    #define a function to rank the molecules by their Tanimoto score against the input molecule.
    def tanimoto_ranker(query): #query must be a SMILES string
        df_mol = pd.read_csv('./TF_DB_clean.csv')
        simil = []
        for i in df_mol['SMILES']:
            simil.append(float(tanimoto_calc(query,i)))
        df_mol['Tanimoto_score_vs_query'] = simil
        df_mol = df_mol.sort_values(by=['Tanimoto_score_vs_query'], ascending=False).reset_index(drop=True)
        return df_mol

    # user input text boxes for SMILES
    user_input_mol = st.text_input("Paste molecule in SMILES format")
    if user_input_mol != "":
        st.write(tanimoto_ranker(user_input_mol))
    else:
        st.write('Please, paste a molecule string in the box above')
