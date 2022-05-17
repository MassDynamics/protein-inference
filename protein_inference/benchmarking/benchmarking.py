import pandas as pd
import numpy as np
import plotly_express as px
import os
import pickle
import matplotlib.pyplot as plt 

from upsetplot import UpSet
import seaborn as sns

from protein_inference.benchmarking.entrapment_benchmark import EntrapmentBenchmark
from protein_inference.benchmarking.benchmark_percolator_inference import ProcessPercolatorInference
from protein_inference.benchmarking.benchmark_percolator_inference import EvaluatePercolatorInference


class Benchmarker():

    def run(self, exp_home,
            true_fastas=[], false_fastas=[], uniprot = False):

        entrapment = true_fastas and false_fastas

        (mq_protein_table,
         mq_peptide_table,
         rep_protein_table,
         rep_peptide_table,
         per_protein_table) = self.collect_tables(exp_home)

        # get entrapment ids
        if entrapment:
            true_targets, true_decoys = self.collect_sequences(
                exp_home, true_fastas, false_fastas)

        # Inference comparison
        mq_inferred = BinaryInferenceComparison().get_maxquant_inferred(mq_protein_table)
        rep_inferred = BinaryInferenceComparison().get_reprisal_inferred(rep_protein_table, uniprot)
        per_inferred = BinaryInferenceComparison().get_percolator_inferred(per_protein_table, uniprot)

        BinaryInferenceComparison().compare_inference(
            [mq_inferred, rep_inferred, per_inferred])

        if entrapment:
            BinaryInferenceComparison().compare_inference_ground_truth([mq_inferred, rep_inferred,
                                                                        per_inferred],
                                                                       true_targets=true_targets,
                                                                       true_decoys=true_decoys)

        # percolator comparison:
        if entrapment:
            EvaluatePercolatorInference().plot_percolator_reprisal_predictions(
                rep_protein_table, per_protein_table, true_targets, true_decoys)
        else:
            EvaluatePercolatorInference().plot_percolator_reprisal_predictions(
                rep_protein_table, per_protein_table)

        GroupingComparison().heatmap_peptide_classifications(mq_peptide_table,
                                                             rep_peptide_table)
        return 1

    def collect_sequences(self, exp_home, true_fastas, false_fastas):

        fasta_home = os.path.join(exp_home, "sequences")
        true_targets = []
        for fasta in true_fastas:
            targets = EntrapmentBenchmark().get_fasta_ids(os.path.join(fasta_home, fasta))
            true_targets = true_targets + targets

        true_decoys = []
        for fasta in false_fastas:
            decoys = EntrapmentBenchmark().get_fasta_ids(os.path.join(fasta_home, fasta))
            true_decoys = true_decoys + decoys

        return true_targets, true_decoys

    def collect_tables(self, exp_home):

        # create paths
        mq_home = os.path.join(exp_home, "MaxQuant")
        rep_home = os.path.join(exp_home, "REPRISAL")
        per_home = os.path.join(exp_home, "MDDiscovery")

        # read files
        mq_protein_table = pd.read_csv(os.path.join(
            mq_home, "reprisal_format", "protein_table.csv"), index_col=False)
        mq_peptide_table = pd.read_csv(os.path.join(
            mq_home, "reprisal_format", "peptide_table.csv"), index_col=False)

        rep_protein_table = pd.read_csv(
            os.path.join(rep_home, "protein_table.csv"))
        rep_peptide_table = pd.read_csv(
            os.path.join(rep_home, "peptide_table.csv"))

        per_protein_table = ProcessPercolatorInference().load_protein_table(
            os.path.join(per_home, "percolator.target.proteins.txt"))

        return mq_protein_table, mq_peptide_table, rep_protein_table, rep_peptide_table, per_protein_table

    def collect_networks(self, exp_home):

        # get mq_networks
        mq_network_path = os.path.join(
            exp_home, "MaxQuant", "reprisal_format", "mq_networks.p")
        mq_networks = pickle.load(open(mq_network_path, "rb"))

        # get rep_networks
        rep_network_path = os.path.join(
            exp_home, "REPRISAL", "target_networks.p")
        rep_networks = pickle.load(open(rep_network_path, "rb"))

        return mq_networks, rep_networks

class BinaryInferenceComparison():

    def get_maxquant_inferred(self, maxquant_PI_df):
        proteins = maxquant_PI_df[maxquant_PI_df["q-value"] <
                                  0.01].protein_id.to_list()
        return proteins

    def get_reprisal_inferred(self, reprisal_PI_df,  uniprot = False):
        if uniprot:
            reprisal_PI_df = self.uniprot_convert_pid(reprisal_PI_df.copy())
        proteins = reprisal_PI_df[reprisal_PI_df["q-value"] <
                                  0.01].protein_id.to_list()
        return proteins

    def get_percolator_inferred(self, percolator_PI_df, uniprot = False):
        if uniprot:
            percolator_PI_df = self.uniprot_convert_pid(percolator_PI_df.copy(),True)
        proteins = percolator_PI_df[percolator_PI_df["q-value"]
                                    < 0.01].ProteinId.to_list()
        return proteins

    def compare_inference(self, list_of_lists, df_names=["MaxQuant", "REPRISAL", "Percolator"], title="Inference Comparison"):

        n = len(list_of_lists)
        results = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                set_1 = set(list_of_lists[i])
                set_2 = set(list_of_lists[j])
                both = len(set_1.intersection(set_2))
                results[i, j] = both

        df = pd.DataFrame(results, columns=df_names, index=df_names)
        fig = px.imshow(df, zmin=0, labels=dict(color="Proteins Inferred"))
        fig.update_layout(title=dict(text=title))
        fig.show()

        return

    def compare_inference_ground_truth(self, list_of_lists, df_names=["MaxQuant", "REPRISAL", "Percolator"], true_targets=[], true_decoys=[]):
        '''lists of proteins coming out of each methodology'''
        true_list_of_lists = [set(list).intersection(
            set(true_targets)) for list in list_of_lists]
        self.compare_inference(true_list_of_lists, title="True Inferences")
        false_list_of_lists = [set(list).intersection(
            set(true_decoys)) for list in list_of_lists]
        self.compare_inference(false_list_of_lists, title="False Inferences")

        return

    def uniprot_convert_pid(self, table, per = False):

        if not per:
            exp = r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]{1,2})"
            table["protein_id"] = table.protein_id.str.extract(exp)
        else:
            exp = r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]{1,2})"
            table["ProteinId"] = table.ProteinId.str.extract(exp)

        return table

    def create_upset_plot(self, mq_protein_table, rep_protein_table, per_protein_table, uniprot = False):
        if uniprot:
            exp = r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]{1,2})"
            rep_protein_table["protein_id"] = rep_protein_table.protein_id.str.extract(exp)
            exp = r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]{1,2})"
            per_protein_table["protein_id"] = per_protein_table.ProteinId.str.extract(exp)
        else:
            per_protein_table["protein_id"] = per_protein_table.ProteinId
        # Inference comparison
        mq_inferred = BinaryInferenceComparison().get_maxquant_inferred(mq_protein_table)
        rep_inferred = BinaryInferenceComparison().get_reprisal_inferred(rep_protein_table)
        per_inferred = BinaryInferenceComparison().get_percolator_inferred(per_protein_table)

        comparison_df = pd.merge(pd.merge(rep_protein_table, per_protein_table, how = "outer"), mq_protein_table, on = "protein_id", how = "outer",suffixes = ("_rep","_mq"))
        comparison_df["mq_inferred"] = comparison_df.protein_id.apply(lambda x: x in mq_inferred)
        comparison_df["rep_inferred"] = comparison_df.protein_id.apply(lambda x: x in rep_inferred)
        comparison_df["per_inferred"] = comparison_df.protein_id.apply(lambda x: x in per_inferred)
        comparison_df = comparison_df[comparison_df.mq_inferred | comparison_df.rep_inferred | comparison_df.per_inferred]
        #comparison_df = comparison_df[["protein_id", "mq_inferred","rep_inferred","per_inferred"]]
        comparison_df = comparison_df.set_index(["mq_inferred","rep_inferred","per_inferred"])

        upset = UpSet(comparison_df, intersection_plot_elements=3,show_percentages=True)
        #upset.add_catplot(value='median_value', kind='strip', color='blue')
        upset.add_catplot(value='total_peptides_rep', kind='box', color='black')
        upset.add_catplot(value='total_peptides_mq', kind='box', color='black')
        fig = plt.figure(figsize=(16,4))
        upset.plot(fig = fig)
        plt.show()

class GroupingComparison():

    def heatmap_peptide_classifications(self, mq_peptide_table, rep_peptide_table):
        '''Look at non duplicate peptides (remove modifications).'''
        df = pd.merge(self.get_peptides_table_from_mod(mq_peptide_table),
                      self.get_peptides_table_from_mod(rep_peptide_table),
                      on="sequence", suffixes=("_mq", "_rep"))
        df = df[["unique_mq", "unique_rep", "razor_mq", "razor_rep"]]
        corr_df = df.corr()
        fig = px.imshow(corr_df,  labels=dict(color="Correlations"))
        fig.update_layout(title=dict(text="Peptide Classification Comparison"))
        fig.show()
        return

    def get_peptides_table_from_mod(self, peptide_table):

        exp = r'n{0,1}\[[\-{0.1}[0-9]+\]'
        peptide_table["sequence"] = peptide_table.sequence_modified.str.replace(
            exp, '')
        peptide_table = peptide_table[["sequence", "razor", "unique", "major"]]
        peptide_table = peptide_table.drop_duplicates()
        peptide_table = peptide_table.sort_values("sequence")

        return peptide_table
