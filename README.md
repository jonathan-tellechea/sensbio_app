
# About Sensbio

Sensbio-app is a simple webserver that allows users to identify possible transcription factors (TFs) inducible by small chemical compounds.

# Required packages

- python 3
- streamlit
- rdkit
- seaborn
- biopython

In addition, a local ncbi-blast+ installation.

# Usage

Using sensbio-app is really simple. The application is separated in three different pages:

- **Data visualization page:** use it to get to know the data in the sensbio dataset.
    - Visualize and export the whole database.
    - Visualize the molecular variety of the database.
    - Visualize the genetic variety of the database.
- **Molecular prediction page:** use to identify possible TFs triggered by your input molecule.
- **TF sequence prediction page:** use to TFs in the database similar to your input protein sequence.

# License

This project is licensed under the terms of the MIT license.
