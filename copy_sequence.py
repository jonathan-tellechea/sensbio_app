import streamlit as st

import pandas as pd
from Bio.Emboss.Applications import NeedleCommandline
import os

def app():
    st.markdown("""
    # Sequence similarity tool

    Use this page to find similar entries in the database using your target protein sequence as input.

    """)
    df_seq = pd.read_csv('./TF_DB_clean.csv')
    seqs = df_seq['AA_sequence']
    with open('database_file.faa','w') as data_file:
        for n,seq in enumerate(seqs):
            data_file.write(f'>{n}\n{seq}\n')

    user_input_seq = st.text_input("Paste sequence string")
    if user_input_seq:
        with open('query.faa','w') as infile:
            infile.write('>query'+'\n'+user_input_seq+'\n')

        needle_cline = NeedleCommandline(asequence="query.faa", bsequence="database_file.faa",gapopen=10, gapextend=0.5,aformat="score", outfile="needle.txt")
        stdout,stderr=needle_cline()

        with open("needle.txt", "r") as infile:
            with open("temp.txt", "w") as output:
                # iterate all lines from file
                for line in infile:
                    # if line starts with substring '#' then don't write it in temp file
                    if not line.strip("\n").startswith('#'):
                        output.write(line)

        # replace file with original name
        os.replace('temp.txt', 'needle.txt')

        ndl_results = pd.read_csv('needle.txt', sep=" ", header=None)
        ndl_results.columns = ["query", "sequence", "alignment", "score"]
        ndl_results['score'] = ndl_results['score'].str.replace('(\(|\))', '').astype(float)
        df_seq['score'] = ndl_results['score']
        df_seq = df_seq.sort_values(by=['score', 'Molecule', 'TF'], ascending=[False, True, True]).reset_index(drop=True) #sort the results by some parameters to not let it to chance
        st.write(df_seq)
    else:
        st.write('Please, paste a sequence string in the box above')
