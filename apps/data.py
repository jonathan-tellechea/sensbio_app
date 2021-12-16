import streamlit as st
import pandas as pd

from PIL import Image
image = Image.open('output.png')

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np
import networkx as nx
import itertools

import plotly.graph_objects as go

def smi_dm(queries):
    lst = []
    n=0
    for smi1, smi2 in itertools.product(queries, repeat = 2):
        dist = 1 - tanimoto_calc(smi1, smi2)
        lst.append(dist)
        n=n+1
        print(str(n)+'/'+str(len(queries)**2), end='\r')
    for i in range(0, len(lst), len(queries)):
        yield lst[i:i + len(queries)]

def tanimoto_calc(smi1, smi2):
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 5, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 5, nBits=2048)
    sim = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
    return sim

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

    st.image(image)

    @st.cache(allow_output_mutation=True)
    def graphCreator(dataframe):
        names = df['Molecule'].drop_duplicates()
        molecules = df['SMILES'].drop_duplicates()
        matrix = np.array([i for i in smi_dm(molecules)])
        dt = [('len', float)]
        A = matrix.view(dt)
        array = A.astype(float)
        threshold = 0.75
        array[array > threshold] = 0.0
        G = nx.from_numpy_matrix(array)
        G = nx.relabel_nodes(G, dict(zip(range(len(G.nodes())), names)))
        return G

    def plotly_graph(graph):
        pos = nx.drawing.nx_agraph.graphviz_layout(graph, prog='neato')
        edge_x = []
        edge_y = []
        for edge in graph.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.append(x0)
            edge_x.append(x1)
            edge_x.append(None)
            edge_y.append(y0)
            edge_y.append(y1)
            edge_y.append(None)

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='#888'),
            hoverinfo='none',
            mode='lines')

        node_x = []
        node_y = []
        for node in graph.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)

        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            hoverinfo='text',
            marker=dict(
                showscale=True,
                colorscale='YlGnBu',
                reversescale=True,
                color=[],
                size=10,
                colorbar=dict(
                    thickness=15,
                    title='Molecule Connections',
                    xanchor='left',
                    titleside='right'
                ),
                line_width=2))
        node_adjacencies = []
        node_text = []
        for node, adjacencies in enumerate(graph.adjacency()):
            node_adjacencies.append(len(adjacencies[1]))
            #node_text.append('# of connections: '+str(len(adjacencies[1])))

        node_trace.marker.color = node_adjacencies
        nodes = [node for node in graph.nodes()]
        #node_trace.text = node_text
        node_trace.text = nodes
        fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                    autosize=True,
                    width=900,
                    height=900,
                    title='<br>Database molecule network',
                    titlefont_size=16,
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    annotations=[ dict(
                        showarrow=False,
                        xref="paper", yref="paper",
                        x=0.005, y=-0.002 ) ],
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
        #fig.show()
        st.plotly_chart(fig, use_container_width=True)

    if st.button('Load molecule network'):
        G = graphCreator(df)
        plot = plotly_graph(G)
        st.set_option('deprecation.showPyplotGlobalUse', False)
        st.pyplot(plot)