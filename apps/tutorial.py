import streamlit as st

def app():
    st.markdown("""
    # Sensbio tutorial

    This page gives a quick overview on how to use the two tools presented in Sensbio.

    ## Data visualization

    This page helps the user visualize and understand the database.
    
    It shows the whole dataset with the following columns:
    - Molecule: common name of the molecule.
    - SMILES: compound in SMILES format
    - InChI: compound in InChI format
    - Species: organism where the TF is found
    - TF: name of the TF gene
    - Bibliographic_ref: PMID or DOI information of the reference describing the TF-molecule interaction.
    - Database_ref: this provides the database reference if the data comes from a database.
    - NCBI_Accession: NCBI acc ID
    - UniProt: UniProt ID
    - AA_sequence: protein sequence of the TF.
    - Pathways: most significative metabolic pathways where the molecule plays a role.

    The user can download the whole dataset as CSV by clicking "Download data as CSV.

    The application generates a barchart describing the pathway distribution of the whole dataset.

    A network graph of the whole database can be visualized by clicking "Load molecule network".

    ## Molecular similarity
    This tool calculates the similarity of the target SMILES with all the molecules present in the database. Apart from the columns already present in the data visualization dataset, the application adds:
    - Tanimoto score vs query: the internal algorithm uses Tanimoto score to rank the dataset by the similarity of the query.
    
    The output shown in the application are the 100 closest entries. The user can download the whole sorted dataset by clicking "Download data as CSV".

    Additionally, the application generates a barchart of the shown dataset describing the pathway distribution of the dataset.


    ## Sequence similarity
    This tool identifies significative BLAST hits in the database using the user's input as query. Apart from the columns already present in the data visualization dataset, the application adds the following information coming from BLAST:
    - id_pc: identity percentage (query/subject)
    - e_value: E-value of the alignment
    - bit_score: alignment score

    As with the other tools, the application generates a barchart of the presented dataset describing the pathway distribution.

    """)