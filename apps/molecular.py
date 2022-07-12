import streamlit as st
import pandas as pd
import numpy as np
import tensorflow as tf

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

import seaborn as sns
import collections
import matplotlib.pyplot as plt

from precalc import smi_cod


def barchartDataframe(dataframe):
    #Graphing the dataset
    molecules = dataframe.drop_duplicates(subset=['Molecule'])
    temp = list(molecules['Pathways'].astype(str))
    temp = [x.split(';') for x in temp if x != 'nan']
    paths = []
    for i in temp:
        for x in i:
            if x != '':
                paths.append(x)
    c = collections.Counter(paths)
    c = c.most_common()
    x = [i[0] for i in c]
    y = [i[1] for i in c]
    p = [i/sum(y)*100 for i in y]
    
    fig, ax = plt.subplots()
    ax = sns.set(rc={'figure.figsize':(60,30)})
    ax = sns.barplot(x=x , y=p)
    ax.axes.set_title("Pathway distribution in the selected dataset",fontsize=120)
    ax.set_ylabel("%",fontsize=120)
    ax.set_xticklabels(ax.get_xticklabels(),rotation = 90)
    ax.tick_params(labelsize=60)
    st.pyplot(fig)

def app():
    st.markdown("""
    # Molecular similarity tool

    Use this page to find similar entries in the database using your target compound in SMILES format as input.

    The final score obtained in the table is the result of adding a 10% of the affinity prediction score to the Tanimoto score.    

    """)
    
    @st.experimental_memo
    def load_files():
        df_mol = pd.read_csv('./TF_DB_clean_pathway.csv')
        aamatrix = np.load('./sequences_matrix.npz')['arr_0']  
        pred_model = tf.keras.models.load_model('./final_model')
        
        return df_mol, aamatrix, pred_model
    
    df_mol, aamatrix, pred_model = load_files()
    
    
    # validate SMILE input
    def chem_validator(query):
        mol = Chem.MolFromSmiles(query)
        if mol is None:
            return False
        return True
    
    #defines a function to calculate the tanimoto similarity between two molecules in SMILES format (from https://medium.com/data-professor/how-to-calculate-molecular-similarity-25d543ea7f40).
    def tanimoto_calc(smi1, smi2):
        mol1 = Chem.MolFromSmiles(smi1)
        mol2 = Chem.MolFromSmiles(smi2)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 5, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 5, nBits=2048)
        sim = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
        return sim
    
    #define a function to predict affinity between a molecule and a sequence
    def model_prediction(aamatrix, fp, model):            
        pred_score = model.predict([aamatrix, fp])            
        return pred_score

    #define a function to rank the molecules by their Tanimoto score against the input molecule.
    def tanimoto_ranker(query, model): #query must be a SMILES string
        simil = []
        for i in df_mol['SMILES']:
            simil.append(float(tanimoto_calc(query,i)))
        df_mol['Tanimoto_score_vs_query'] = simil
        
        fp = smi_cod(query)
        fp = np.repeat(fp, np.shape(aamatrix)[0], axis=0)
        
        preds = model_prediction(aamatrix,fp,model)
        affin = preds[:,0].tolist()
        
        #final score
        score = (np.array(simil)+0.1*np.array(affin)).tolist()
            
        df_mol['Affinity prediction'] = affin
        df_mol['Final score'] = score 
        
        df_mol_sorted = df_mol.sort_values(by=['Tanimoto_score_vs_query'], ascending=False).reset_index(drop=True)
        return df_mol_sorted
        

    # user input text boxes for SMILES
    user_input_mol = st.text_input("Paste molecule in SMILES format")
    if user_input_mol != "" and chem_validator(user_input_mol):
        df = tanimoto_ranker(user_input_mol, pred_model)
        # df_head = df.head(100)
        st.write(df)

        @st.cache
        def convert_df(dataframe):
            # IMPORTANT: Cache the conversion to prevent computation on every rerun
            return dataframe.to_csv().encode('utf-8')

        csv = convert_df(df)

        st.download_button(
            label="Download data as CSV",
            data=csv,
            file_name='molecular_prediction_results.csv',
            mime='text/csv'
        )
        # barchartDataframe(df_head)
        barchartDataframe(df)
    else:
        st.write("Please, paste your target molecule in SMILES format in the box above or paste some example molecules:")
        st.code('''Examples:
        naringenin: O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c21
        vanillate: COC1=C(C=CC(=C1)C(=O)O)O
        uracil: Oc1ccnc(O)n1''')
    
    
