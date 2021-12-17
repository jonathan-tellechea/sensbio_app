import streamlit as st
from multiapp import MultiApp
from apps import data, molecular, sequence, tutorial# import your app modules here


import pandas as pd

app = MultiApp()

st.markdown("""
# Sensbio: An online server for biosensor design

The application is divided into three interactive pages.

- **Data visualization:** designed to help to better understand the sensbio tool.
- **Molecular prediction:** tool: the user can find molecules similar to their input and find candidate TF where to start the design of new biosensor circuits specific for the target molecule.
- **Sequence prediction tool:** the user can find similar protein sequences similar to the ones already described in the database and the ligand molecule that triggers them.

""")

# Add all application here
app.add_app("Data visualization", data.app)
app.add_app("Molecular similarity", molecular.app)
app.add_app("Sequence similarity", sequence.app)
app.add_app("Tutorial", tutorial.app)
# The main app
app.run()

st.markdown("""**Attribution note:** The multi-page functionality is using the [streamlit-multiapps](https://github.com/upraneelnihar/streamlit-multiapps) framework developed by [Praneel Nihar](https://medium.com/@u.praneel.nihar).
""")
