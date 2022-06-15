import streamlit as st

from Bio.Blast.Applications import NcbiblastpCommandline
import pandas as pd
import numpy as np
import tensorflow as tf

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

import seaborn as sns
import collections
import matplotlib.pyplot as plt

from precalc import seq_cod

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
    # Sequence similarity tool
    Use this page to find similar entries in the database using your target protein sequence as input.
    """)
    
    database = pd.read_csv('./TF_DB_clean_pathway.csv')
    fpmatrix = np.load('./fingerprints_matrix.npy')
    pred_model = tf.keras.models.load_model('./final_model')
    
    #define a function to predict affinity between a molecule and a sequence
    #define a function to predict affinity between a molecule and a sequence
    def model_prediction(aa, fpmatrix, model):            
        pred_score = model.predict([aa, fpmatrix])            
        return pred_score
    
    
    # st.write("Use these options to decide if predictions are performed over all" 
    #          " the database. Otherwise, they will be performed only on the"
    #          " first 100/500/1000.")
    
    # check = st.radio("Select an option:",
    #                  ("Perform predictions over all the database",
    #                   "Perform predictions over the first 100",
    #                   "Perform predictions over the first 500",
    #                   "Perform predictions over the first 1000"))
    
    # if check == "Perform predictions over all the database":
    #     flag = True
    # elif check == "Perform predictions over the first 100":
    #     flag = False
    #     num_check = 100
    # elif check == "Perform predictions over the first 500":
    #     flag = False
    #     num_check = 500
    # elif check == "Perform predictions over the first 1000":
    #     flag = False
    #     num_check = 1000
    
        
    user_input_seq = st.text_input("Paste sequence string")
    validate = False
    if user_input_seq:
        try:
            blastp_cline = NcbiblastpCommandline(db="BLASTdb", outfmt="6 sseqid pident evalue bitscore") #Blast command
            out, err = blastp_cline(stdin=user_input_seq)
            validate = True
        except:
            print("Wrong sequence format!", user_input_seq)
            validate = False
    if user_input_seq and validate:
        results = [i.split('\t') for i in out.splitlines()]
        for i in results:
            if '|' in i[0]:
                i[0] = i[0].split('|')[1]

        blast_df = pd.DataFrame(data=results, columns=['NCBI_Accession', 'id_pc','e_value','bit_score'])
        #blast_df['NCBI_Accession'] = blast_df['NCBI_Accession'].apply(lambda x : x.split('|')[1])
        blast_df['bit_score'] = pd.to_numeric(blast_df['bit_score'])

        df = database.merge(blast_df).sort_values(by=['bit_score'], ascending=False).reset_index(drop=True)
        
        aa = seq_cod(user_input_seq)
        aa = np.repeat(aa, np.shape(fpmatrix)[0], axis=0)
        
        preds = model_prediction(aa,fpmatrix,pred_model)
        affin = preds[:,0].tolist()
            
        df['Affinity prediction'] = affin
        
        st.write(df)

        @st.cache
        def convert_df(dataframe):
            # IMPORTANT: Cache the conversion to prevent computation on every rerun
            return dataframe.to_csv().encode('utf-8')

        csv = convert_df(df)

        st.download_button(
            label="Download data as CSV",
            data=csv,
            file_name='sequence_prediction_results.csv',
            mime='text/csv'
        )
        barchartDataframe(df)
    else:
        st.write("Please, paste your target protein sequence in the box above or paste some example sequences:")
        st.text_area('LacI (E. coli)', 'MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQSLLIGVATSSLALHAPSQIVAAIKSRADQLGASVVVSMVERSGVEACKAAVHNLLAQRVSGLIINYPLDDQDAIAVEAACTNVPALFLDVSDQTPINSIIFSHEDGTRLGVEHLVALGHQQIALLAGPLSSVSARLRLAGWHKYLTRNQIQPIAEREGDWSAMSGFQQTMQMLNEGIVPTAMLVANDQMALGAMRAITESGLRVGADISVVGYDDTEDSSCYIPPLTTIKQDFRLLGQTSVDRLLQLSQGQAVKGNQLLPVSLVKRKTTLAPNTQTASPRALADSLMQLARQVSRLESGQ')
        st.text_area('AseR (B. subtilis)', 'MTIDVAAMTRCLKTLSDQTRLIMMRLFLEQEYCVCQLVDMFEMSQPAISQHLRKLKNAGFVNEDRRGQWRYYSINGSCPEFDTLQLILHQIDQEDELLNHIKQKKTQACCQ')
