import streamlit as st
import pandas as pd

from PIL import Image
dist = Image.open('output.png')
netw = Image.open('network.png')

def app():

    st.markdown("""
    # Data visualization

    Use this page to visualize the database information and understand what is in it.
    """)
    

    data_load_state = st.text('Loading data...')

    # Load the dataframe.
    df = pd.read_csv('./TF_DB_clean_pathway.csv')
    st.write(df)
    # Notify the reader that the data was successfully loaded.
    data_load_state.text('Loading data...done! This is the whole dataset')

    @st.cache
    def convert_df(dataframe):
        return dataframe.to_csv().encode('utf-8')
    
    csv = convert_df(df)
    
    st.download_button(
        label="Download data as CSV",
        data=csv,
        file_name='database.csv',
        mime='text/csv'
    )

    st.image(dist)
    st.image(netw)
