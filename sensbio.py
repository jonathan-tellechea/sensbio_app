import streamlit as st
from multiapp import MultiApp
from apps import data, molecular, sequence # import your app modules here

app = MultiApp()

st.markdown("""
# Sensbio web application

The application is divided in three interactive pages.

- Data visualization:
    - Designed to help to better understand the sensbio database.
- Molecular prediction tool:
    - The user can find molecules similar to their input and find candidate TF where to start the design of new biosensor circuits specific for the target molecule.
- Sequence prediction tool
    - The user can find similar protein sequence to the ones already described in the database and the ligand molecule that triggers them.

**Attribution note:** The multi-page functionality is using the [streamlit-multiapps](https://github.com/upraneelnihar/streamlit-multiapps) framework developed by [Praneel Nihar](https://medium.com/@u.praneel.nihar).

""")

# Add all your application here
app.add_app("Data visualization", data.app)
app.add_app("Molecular similarity", molecular.app)
app.add_app("Sequence similarity", sequence.app)
# The main app
app.run()
