import streamlit as st
import pandas as pd

def app():

    st.markdown("""
    # Data visualization

    Use this page to visualize the database information and understand what is in it.

    """)
    data_load_state = st.text('Loading data...')

    # Load the dataframe.
    df = pd.read_csv('./TF_DB_clean.csv')
    st.write(df)
    # Notify the reader that the data was successfully loaded.
    data_load_state.text('Loading data...done! This is the whole dataset')
