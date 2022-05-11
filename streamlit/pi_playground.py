import streamlit as st

import numpy as np
import pandas as pd
import pickle
import argparse
from PIL import Image 

# imports for PI
import sys
import os
from pandas import read_csv
import plotly.figure_factory as ff
import plotly.express as px
from io import BytesIO

sys.path.insert(1,"..")
sys.path.insert(1,"../protein_inference")

from protein_inference.table_maker import TableMaker
from protein_inference.network_grapher import NetworkGrapher

# pass some args
def parse_args(args):
    parser = argparse.ArgumentParser('Data Diagnostics')
    parser.add_argument('-f', '--folder', dest = "pi_folder", 
            help='Folder containing protein inference output', required=False)
    return parser.parse_args(args)

args = parse_args(sys.argv[1:])

# app
st.set_page_config(layout="wide")

with st.sidebar:
    background = Image.open("protein_inference_logo_v2.png")
    st.image(background, use_column_width=True)
    st.title("Protein Inference PlayGround")
    st.subheader("Joseph Bloom - Mass Dynamics 2021")

    '''
    Welcome to Protein Inference Playground!
    
    This tool is designed to help you explore protein inference results produced with the PI python package. 
    Protein Inference is the process of inferring proteins present in a sample based on the peptide spectral matches found in a Mass Spectrometry experiment.

    Use the Navigation at the top of the page to choose between:
    - Experiment Summary (This page has your protein inference results) and some summary statistics. 
    - Quality Control (This page has some diagnostic graphs which may help you understand your results.)
    - Network Visualization (This page provides force-directed network graphics to visualize results.)
    '''

    '''
    ---


    Did you enjoy using the protein inference playground? 

    [Let us know on twitter!](https://twitter.com/massdynamicsco)

    Want to give us feedback? Reach out to Mass Dynamics here [here](https://www.massdynamics.com/get-in-touch) or to the authors directly at [here](joseph@massdynamics.com)


    '''

    st.image("new_md_logo.png")
# import libraries
import streamlit as st

@st.cache(allow_output_mutation=True)
def load_protein_inference_data(path):
    with st.spinner(text = "Loading Protein Inference Results"):
        target_protein_table = pd.read_csv(os.path.join(path,"reprisal.target.proteins.csv")).drop(["ProteinGroupId", "FDR"], axis = 1)
        target_peptide_table = pd.read_csv(os.path.join(path,"reprisal.target.peptides.csv"))
        decoy_protein_table = pd.read_csv(os.path.join(path,"reprisal.decoy.proteins.csv"))
        decoy_peptide_table = pd.read_csv(os.path.join(path,"reprisal.decoy.peptides.csv"))
        target_networks = pickle.load(open(os.path.join(path,"target_networks.p"), "rb"))
    return target_protein_table, target_peptide_table, decoy_protein_table, decoy_peptide_table, target_networks

#    st.info("No protein inference output folder has been defined.")
#    st.info("You can define it by calling: streamlit run pi_playground.py -- --folder /path/to/your/pi/output/")

if args.pi_folder is not None:
    target_protein_table, target_peptide_table, \
            decoy_protein_table, decoy_peptide_table, \
            target_networks = load_protein_inference_data(args.pi_folder)
    min_tda_score = target_protein_table[target_protein_table["q-value"] < 0.01].score.min()
    navigation = st.radio('Page Selection', ["Experiment Summary", "Quality Control", "Network Visualization"], index = 1)
else: 
    navigation = "experiment_upload"



if navigation == "experiment_upload":
    all_files_uploaded = False
    needed_files = {"reprisal.target.proteins.csv", "reprisal.target.peptides.csv", "reprisal.decoy.proteins.csv", "reprisal.decoy.peptides.csv", "target_networks.p"}

    st.write("Please upload all the files.")
    uploaded_files = st.file_uploader("Please upload the protein and peptide, target and decoy tables here:", accept_multiple_files=True, key = 102)
    for uploaded_file in uploaded_files:     
        
        if uploaded_file.name == "reprisal.target.proteins.csv":
            uploaded_table = pd.read_csv(uploaded_file)
            target_protein_table = uploaded_table.drop(["ProteinGroupId", "FDR"], axis = 1)
            needed_files = needed_files - {"reprisal.target.proteins.csv"}
        elif uploaded_file.name == "reprisal.target.peptides.csv":
            uploaded_table = pd.read_csv(uploaded_file)
            target_peptide_table = uploaded_table
            needed_files = needed_files - {"reprisal.target.peptides.csv"}
        elif uploaded_file.name == "reprisal.decoy.proteins.csv":
            uploaded_table = pd.read_csv(uploaded_file)
            decoy_protein_table = uploaded_table
            needed_files = needed_files - {"reprisal.decoy.proteins.csv"}
        elif uploaded_file.name =="reprisal.decoy.peptides.csv":
            uploaded_table = pd.read_csv(uploaded_file)
            decoy_peptide_table = uploaded_table
            needed_files = needed_files - {"reprisal.decoy.peptides.csv"}
        elif uploaded_file.name == "target_networks.p":
            with open("target_networks.p", "wb") as f:
                f.write(uploaded_file.getbuffer())
            target_networks = pickle.load(open("target_networks.p", "rb"))
            needed_files = needed_files - {"target_networks.p"}
    
    if len(needed_files) == 0:
        all_files_uploaded = True
        st.write("All files uploaded!")

        min_tda_score = target_protein_table[target_protein_table["q-value"] < 0.01].score.min()

        st.write("You can now navigate to the other pages.")
        navigation = st.radio('Page Selection', ["Experiment Summary", "Quality Control", "Network Visualization"], index = 1)

if navigation == "Experiment Summary":
    if (1): # PLACEHOLDER FOR HAVING FINISHED PROCESSING
     
        st.header("Output Tables:")

        searchbox = st.text_input("Search for proteins here. Use ProteinId strings to search. No fancy syntax like regex or anything. Exact or contains matches only.", "")

        st.subheader("Protein Table")
        if searchbox:
            st.write((target_protein_table[
                target_protein_table.ProteinId.str.contains(searchbox) | 
                target_protein_table.indistinguishable.str.contains(searchbox) |
                target_protein_table.subset.str.contains(searchbox)
            
            ].sort_values("score", ascending = False)))
        else:
            st.write(target_protein_table.sort_values("score", ascending = False).head(5))

        st.subheader("Peptide Table")
        st.write(target_peptide_table)

if navigation == "Network Visualization":
    st.header("Network Visualizer:")
    st.subheader("")
    #Visualize networks:
    if "target_protein_table" in locals():
        '''
        Protein Inference algorithms can be interpreted using network diagrams. 
    
        Please note the following:
            - Large nodes are proteins. Small nodes are peptides (or PSMs).
            - Edges are drawn when a protein contains a peptides (may have emitted it).
            - Coloring options include by annotation, 
        
        In the default visualization mode, proteins and peptides are coloured by status (annotation), group, or score.

        For more info please see the current github repository or reach out to the authors. 
    
        '''
        molecule = st.text_input("Which molecule do you want to visualize? (copy and paste protein id's from the protein table above)", 
                                    target_protein_table.sort_values("non_unique").ProteinId.to_list()[-1])
        pn = TableMaker().find_molecule(target_networks, molecule) 

        group_option = st.selectbox("Select a color schema for the force-directed problem network visualization. Selecting 'all' will create all options below.", ["all", "status", "group", "score"])

        if group_option != "all":
            NetworkGrapher().draw(pn, group_option, size = [800,800])
            st.components.v1.html(open("nx.html", "r").read(), width=800, height=800)
        
        else:
            width = 400
            height = 400
            col1, col2, col3 = st.columns(3)
            NetworkGrapher().draw(pn, "status", name = "nx", size = [width,height])
            NetworkGrapher().draw(pn, "group", name = "nx1", size = [width,height])
            NetworkGrapher().draw(pn, "score", name = "nx2", size = [width,height])
            with col1:
                st.components.v1.html(open("nx.html", "r").read(), width=width, height=height)
                st.download_button('Download', open("nx.html", "r").read(), file_name = "problem_network.html")
            with col2:
                st.components.v1.html(open("nx1.html", "r").read(), width=width, height=height)
                st.download_button('Download', open("nx1.html", "r").read(), file_name = "problem_network.html")
            with col3:
                st.components.v1.html(open("nx2.html", "r").read(), width=width, height=height)
                st.download_button('Download', open("nx2.html", "r").read(), file_name = "problem_network.html")
        
        protein_table_option = st.checkbox("Show Protein Table", key = 1)
        peptide_table_option = st.checkbox("Show Peptide Table", key = 2)

        if protein_table_option:
            st.subheader("Protein Table")
            tmp = TableMaker().get_protein_table(pn)
            tmp["Inferred"] = tmp["score"].apply(lambda x: x > min_tda_score)
            st.write(tmp)

        if peptide_table_option:
            st.subheader("Peptide Table")
            st.write(TableMaker().get_peptide_table(pn))
    else: 
        st.info("I'm thinking about my weekend plans while waiting...")

if navigation == "Quality Control":
    st.subheader("QC Figures")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Decoy - Target Distribution of Log Scores")
        '''
        The following overlayed histograms demonstrate the seperation of scores corresponding to known decoys and target proteins. 
        
        Many protein inference strategies, including REPRISAL, with calculate the FDR of a protein score as the 
        ratio of decoys which attain an equal to or greater score to targets which attain an equal or greater score. 

        '''
        if "target_protein_table" in locals():
            
            decoy_scores = decoy_protein_table.score.apply(np.log).to_list()
            target_scores = target_protein_table.score.apply(np.log).to_list()
            all_scores = target_scores + decoy_scores
            label = ["target"] * len(target_scores) + ["decoy"] * len(decoy_scores)
            tmp = pd.DataFrame({"score": all_scores, "label": label})

            fig = px.histogram(
                tmp, x="score", color="label", 
                nbins=50, barmode="overlay",
                labels={"score": "log(score)", "label": "TDA Group"},
                template = "plotly_dark"
            )

            min_tda_score = target_protein_table[target_protein_table["q-value"] < 0.01].score.min()
            fig.add_vline(x=np.log(min_tda_score), line_dash="dashdot", annotation={"text": "TDA Inference threshold ({})".format(round(min_tda_score,3))})
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.subheader("ECDF of Log Scores by Number of Peptides in the Group")
        '''
        The following ECDFs demonstrate the distribution of scores by number of peptides assigned to a protein (among major proteins.)

        This is useful for getinng a sense for whether homology is a significant challenge within this dataset.
        If proteins similtaneously have low scores and many peptides, then protein inference is "harder".

        '''

        tmp = target_protein_table.copy()
        tmp["total_peptides"] = tmp.total_peptides.apply(lambda x: str(x) if x < 5 else "5+")
        tmp["score"] = tmp["score"].apply(lambda x: np.log(x+0.001))

        min_tda_score = target_protein_table[target_protein_table["q-value"] < 0.01].score.min()
        #print(min_tda_score)

        fig = px.ecdf(tmp.sort_values("total_peptides"), x="score", color = "total_peptides", template = "plotly_dark", 
            labels={"total_peptides": "Number of Peptides Total", "score":"log(score)"}, color_discrete_map={"5+": "red", "1": "blue", "2": "green", "3": "orange", "4": "purple", "5": "black"})
        fig.add_vline(x=np.log(min_tda_score), line_dash="dashdot", annotation={"text": "TDA Inference threshold ({})".format(round(min_tda_score,3))})
        st.plotly_chart(fig, use_container_width=True)
