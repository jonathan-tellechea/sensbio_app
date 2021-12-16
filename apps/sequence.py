import streamlit as st

from Bio.Blast.Applications import NcbiblastpCommandline
import pandas as pd

import seaborn as sns
import collections
import matplotlib.pyplot as plt

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
    
    user_input_seq = st.text_input("Paste sequence string")
    if user_input_seq:
        blastp_cline = NcbiblastpCommandline(db="BLASTdb", outfmt="6 sseqid pident evalue bitscore") #Blast command
        out, err = blastp_cline(stdin=user_input_seq)

        results = [i.split('\t') for i in out.splitlines()]
        for i in results:
            if '|' in i[0]:
                i[0] = i[0].split('|')[1]

        blast_df = pd.DataFrame(data=results, columns=['NCBI_Accession', 'id_pc','e_value','bit_score'])
        #blast_df['NCBI_Accession'] = blast_df['NCBI_Accession'].apply(lambda x : x.split('|')[1])
        blast_df['bit_score'] = pd.to_numeric(blast_df['bit_score'])

        df = database.merge(blast_df).sort_values(by=['bit_score'], ascending=False).reset_index(drop=True)
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
