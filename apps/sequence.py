import streamlit as st

from Bio.Blast.Applications import NcbiblastpCommandline
import pandas as pd

def app():
    st.markdown("""
    # Sequence similarity tool

    Use this page to find similar entries in the database using your target protein sequence as input.

    """)
    database = pd.read_csv('./TF_DB_clean.csv')
    
    user_input_seq = st.text_input("Paste sequence string")
    if user_input_seq:
        blastp_cline = NcbiblastpCommandline(db="BLASTdb", outfmt="6 sseqid pident evalue bitscore") #Blast command
        out, err = blastp_cline(stdin=user_input_seq)

        results = [i.split('\t') for i in out.splitlines()]

        blast_df = pd.DataFrame(data=results, columns=['NCBI_Accession', 'id_pc','e_value','bit_score'])
        # blast_df['NCBI_Accession'] = blast_df['NCBI_Accession'].apply(lambda x : x.split('|')[1])
        blast_df['bit_score'] = pd.to_numeric(blast_df['bit_score'])

        scored_db = database.merge(blast_df).sort_values(by=['bit_score'], ascending=False).reset_index(drop=True)
        scored_db
        st.write(scored_db)
    else:
        st.write('Please, paste a sequence string in the box above')
