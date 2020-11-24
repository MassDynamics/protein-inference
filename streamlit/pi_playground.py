import streamlit as st

import numpy as np
import pandas as pd
import pickle

# imports for PI
import sys
import os
from pandas import read_csv
sys.path.insert(1,"..")
sys.path.insert(1,"../protein_inference")

from protein_inference.table_maker import TableMaker
from protein_inference.network_grapher import NetworkGrapher

# app
st.image("thumbnail_graph.png")
st.title("Protein Inference PlayGround")
st.subheader("Joseph Bloom - Mass Dynamics 2020")

'''
Welcome to Protein Inference Playground. This tool is designed to help you understand protein inference results.

Protein Inference is the process of inferring proteins present in a sample based on the peptide spectral matches found in a Mass Spectrometry experiment.

The protein inference workflow used here is composed of the following steps:
   - PSM Filtering by q-value
   - PSM Network Representation
   - Protein/Peptide Grouping and Scoring
   - Estimation of false discovery rates
   - Table and Output Graph Writing

The only inference method demonstrated in this app is MD's own **REPRISAL** (**R**ecursive **Pr**otein **I**nference **S**coring **Al**gorithm).

The data below corresponds to the iPRG2015 dataset which is not designed to test protein inference algorithms with respect to high homology, but 
nevertheless contains many examples of PSM networks that can be visualized. 
'''

@st.cache(allow_output_mutation=True)
def get_output_data(path):
    with st.spinner(text = "Loading Protein Inference Results"):
        target_protein_table = pd.read_csv(os.path.join(path,"reprisal.target.proteins.csv"))
        target_peptide_table = pd.read_csv(os.path.join(path,"reprisal.target.peptides.csv"))
        decoy_protein_table = pd.read_csv(os.path.join(path,"reprisal.decoy.proteins.csv"))
        decoy_peptide_table = pd.read_csv(os.path.join(path,"reprisal.decoy.peptides.csv"))
        target_networks = pickle.load(open(os.path.join(path,"target_networks.p"), "rb"))
    return target_protein_table, target_peptide_table, decoy_protein_table, decoy_peptide_table, target_networks

target_protein_table, target_peptide_table, decoy_protein_table, decoy_peptide_table, target_networks = get_output_data("../example_data/IPRG2015/")

st.header("Output Tables:")
if st.checkbox("Show Protein Table"):
    st.subheader("Protein Table")
    st.write(target_protein_table.sort_values("score", ascending = False))

if st.checkbox('Show Peptide Table'):
    st.subheader("Peptide Table")
    st.write(target_peptide_table)

import plotly.figure_factory as ff
st.header("Decoy - Target Distribution of Log Scoress:")
'''
The following overlayed histograms demonstrate the seperation of scores corresponding to known decoys and target proteins. Target proteins are not necessarily present but decoy proteins are certainly absent, hence some overlap. 

'''
if "target_protein_table" in locals():
    # Add histogram data
    target = np.log(target_protein_table.score+0.01)
    decoy = np.log(decoy_protein_table.score+0.01)

    # Group data together
    hist_data = [target, decoy]
    group_labels = ['Target', 'Decoy']
    fig = ff.create_distplot(hist_data, group_labels, bin_size = [0.3,0.3])
    st.plotly_chart(fig, use_container_width=True)
else:
    st.info("I'm just bored because you haven't uploaded any files yet...")

st.header("Network Visualizer:")
st.subheader("")
#Visualize networks:
if "target_protein_table" in locals():
    '''
    Protein Inference algorithms can be interpreted using network diagrams. 
    
    These networks are composed of protein nodes (larger) and peptide nodes (smaller). Edges are drawn between protein and peptide nodes where the peptide is contained within the protein sequence and we believe the peptide is present in your sample (a Peptide Spectral Match or PSM). 

    In the default visualization mode, proteins and peptides are coloured by their status. 
    
    In the "Colour by Group" mode, proteins and peptides are coloured by their associated group. For a peptide, a group indicated the major protein which we think emitted it. For a protein, a group may mean different things depending on the algorithm that allocated protein groupings. 

    In REPRISAL, protein groups correspond to either the protein itself it is major, or the protein group that explains the last remaining unexplained peptides during iterative solving. In most cases this collapses to being a protein which explains all peptides associated with that protein.

    '''
    molecule = st.text_input("Which molecule do you want to visualize? (copy and paste protein id's from the protein table above)", "sp|P33302|PDR5_YEAST")
    pn = TableMaker().find_molecule(target_networks, molecule) 

    group_option = st.checkbox("Colour by Group", key = 0)
    protein_table_option = st.checkbox("Show Protein Table", key = 1)
    peptide_table_option = st.checkbox("Show Peptide Table", key = 2)

    if group_option:
        NetworkGrapher().draw(pn, "group", size = [800,800])
    else: 
        NetworkGrapher().draw(pn, size = [800,800])

    myfile = open("nx.html", "r")
    data= myfile.read()
    st.components.v1.html(data, width=800, height=800)

    if protein_table_option:
        st.subheader("Protein Table")
        st.write(TableMaker().get_protein_table(pn))

    if peptide_table_option:
        st.subheader("Peptide Table")
        st.write(TableMaker().get_peptide_table(pn))
else: 
    st.info("I'm thinking about my weekend plans while waiting...")


'''
---


Did you enjoy using the protein inference playground? 

[Let us know on twitter!](https://twitter.com/massdynamicsco)

Want to give us feedback? Reach out to us [here](https://www.massdynamics.com/get-in-touch):


'''

st.image("new_md_logo.png")