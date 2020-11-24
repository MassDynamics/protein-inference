## sidebar

st.sidebar.write("Please note that FDR cut off's and alternative solutions are not currently implemented.")
add_selectbox = st.sidebar.selectbox(
    'Which solution Method would you like to use?',
    ('REPRISAL', 'Percolator', 'FIDO')
)

psm_fdr = st.sidebar.multiselect(
    'Select a PSM FDR Cutoff',
    [0.01, 0.05]
)

protein_fdr = st.sidebar.multiselect(
    'Select a Protein FDR Cutoff',
    [0.01, 0.05]
)


@st.cache(allow_output_mutation=True)
def get_data(psms, decoy=0, scoring_method = GreedyAlgorithm):
    psms = PSMsPreprocessor(psms, decoy=decoy).get_processed_psms()
    network = PSMsNetworkGenerator(psms).generate_network()
    problem_networks = PSMsNetworkSplitter(network).split_networks()
    tagged_networks = ProteinInferenceRunner().parallel_apply(
        problem_networks, UniquenessTagger().run)
    solved_networks = ProteinInferenceRunner().parallel_apply(
        tagged_networks, scoring_method().run)
    merged_networks = ProteinInferenceRunner().parallel_apply(
        solved_networks, ProteinMerger().run)
    protein_table = TableMaker().get_system_protein_table(merged_networks)
    peptide_table = TableMaker().get_system_peptide_table(merged_networks)
    return merged_networks, protein_table, peptide_table

    


#File Uploading ( unfinished)
#st.header("Upload your PSM Files:")
##remove warning
#st.set_option('deprecation.showfileUploaderEncoding', False)
#uploaded_file = st.file_uploader("Choose a Target PSMs List.", type=["csv","txt","tsv"], key = 11)
#if uploaded_file is not None:
#    target_df = pd.read_csv(uploaded_file, sep="\t", engine="python")
#uploaded_file = st.file_uploader("Choose a Decoy PSMs List.", type=["csv","txt","tsv"],key = 10)
#if uploaded_file is not None:
#    decoy_df = pd.read_csv(uploaded_file, sep="\t", engine="python")


#try:
#    with st.spinner(text = "Running Protein Inferece..."):
#        tagged_networks, solved_networks, merged_networks, protein_table, peptide_table = get_data(target_df)
#        decoy_networks, decoy_networks, decoy_networks, decoy_protein_table, decoy_peptide_table = get_data(decoy_df,1)
#        target_fdr_table = FalseDiscoveryRateCalculator().tag_FDR(protein_table,
#                                                                decoy_protein_table)
#        st.success("Your data is ready!")
#except:
#        st.warning("No data loaded yet.")


# do processing (initial)
#target_path = "../example_data/IPRG2015/percolator.target.psms.txt"
#decoy_path = "../example_data/IPRG2015/percolator.decoy.psms.txt"

#target_df = pd.read_csv(target_path, sep = "\t")
#decoy_df = pd.read_csv(decoy_path, sep = "\t")

#with st.spinner(text = "Running Protein Inference..."):
#    merged_networks, protein_table, peptide_table = get_data(target_df)
#    _, decoy_protein_table, _ = get_data(decoy_df,1)
#    target_fdr_table = FalseDiscoveryRateCalculator().tag_q_value(protein_table, decoy_protein_table)