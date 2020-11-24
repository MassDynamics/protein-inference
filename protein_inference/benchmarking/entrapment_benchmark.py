from protein_inference.inference import FalseDiscoveryRateCalculator
from protein_inference.protein_inference_runner import ProteinInferenceRunner
import numpy as np
import plotly_express as px
import pandas as pd
import os

class EntrapmentBenchmark():

    def run(self, experiment_home,  positives, negatives):
        true_positives, true_negatives, target_psms, decoy_psms = self.load_data_for_entrapment_benchmarking(
            experiment_home,  positives, negatives)

        target_protein_table = self.benchmark_FDRs(target_psms, decoy_psms, true_negatives)
        self.boxplot_FDR_dif(target_protein_table)
        self.plot_ent_fdr_with_decoy_fdr(target_protein_table)
        self.plot_pos_with_fdr(target_protein_table, true_positives)

        benchmark_df = target_protein_table
        preds = benchmark_df[benchmark_df["q-value"] < 0.01]
        true_positive_preds = preds[preds.ProteinId.apply(
            lambda x: x in true_positives)]
        false_positive_preds =  preds[preds.ProteinId.apply(
            lambda x: x in true_negatives)]

        print("number of true positives: ", len(true_positives))
        print("number of true negatives: ", len(true_negatives))

        print("At 0.01% FDR")
        print("Number of True positives: ", len(true_positive_preds))
        print("Recall: ", len(true_positive_preds)/len(true_positives))
        print("Number of False positives: ", len(false_positive_preds))
        print("FDR: ", len(false_positive_preds)/preds.shape[0])
        
        return target_protein_table

    def benchmark_FDRs(self, known_psms, generated_decoy_psms, true_negatives):

        #known_psms = known_psms[known_psms["percolator q-value"] < 0.001]

        _, true_decoy_psms = self.get_true_and_entrapment_sets(
            known_psms, [], true_negatives)

        target_protein_table, _, _ = ProteinInferenceRunner().get_output(known_psms)
        entrapment_protein_table, _, _ = ProteinInferenceRunner().get_output(true_decoy_psms, 1)
        decoy_protein_table, _, _ = ProteinInferenceRunner().get_output(generated_decoy_psms, 1)

        # tag FDR's.
        target_protein_table = FalseDiscoveryRateCalculator().tag_q_value(target_protein_table,  decoy_protein_table)
        target_protein_table = FalseDiscoveryRateCalculator().tag_q_value(target_protein_table,  entrapment_protein_table, entrapment=True)
       
        target_protein_table["FDR_dif"] = target_protein_table["q-value"] - \
            target_protein_table["q-value-entrapment"]


        return target_protein_table

    def boxplot_FDR_dif(self, benchmark_df):
        benchmark_df = benchmark_df[benchmark_df.entrapmentFDR > 0.15]
        fig = px.box(benchmark_df, y="FDR_dif")
        fig.update_yaxes(range=[-0.2, 0.2])
        fig.show()
        return

    def plot_ent_fdr_with_decoy_fdr(self, benchmark_df):
        benchmark_df = benchmark_df.rename(
            {"FDR": "Decoy FDR", "entrapmentFDR": "Entrapment FDR"}, axis="columns")
        fig = px.scatter(benchmark_df, x="Decoy FDR", y="Entrapment FDR",
                         trendline="ols", template="simple_white")
        fig.show()
        return

    def plot_pos_with_fdr(self, benchmark_df, true_positives):
        true_positive_preds = benchmark_df[benchmark_df.ProteinId.apply(
            lambda x: x in true_positives)]
        fdr_space = np.linspace(0, 0.1, 3000)
        pos_count = [sum(true_positive_preds["q-value"] < fdr) for fdr in fdr_space]
        df = pd.DataFrame(
            {"FDR": fdr_space, "Number of Protein Groups": pos_count})
        fig = px.line(df, x="FDR", y="Number of Protein Groups",
                      template="simple_white", )
        fig.show()
        return

    def get_fasta_ids(self, fasta_file):

        if (type(fasta_file) is list):
            for file in fasta_file:
                positives = []
                positives_file = open(file, 'r')
                for line in positives_file.readlines():
                    if line.startswith(">"):
                        positives.append(line[1:-1])

            return positives
        
        else:
            positives = []
            positives_file = open(fasta_file, 'r')
            for line in positives_file.readlines():
                if line.startswith(">"):
                    positives.append(line[1:-1])

        return positives
    
    def get_true_and_entrapment_sets(self, psms, true_positives, true_negatives):
        true_target_psms = psms[psms["protein id"].apply(
            lambda x: x in true_positives)]
        true_decoy_psms = psms[psms["protein id"].apply(
            lambda x: x in true_negatives)]
        return true_target_psms, true_decoy_psms

    def load_data_for_entrapment_benchmarking(self, experiment_home, positives, negatives):

        if type(negatives) is list:
            true_negatives = self.get_fasta_ids(
                [os.path.join(experiment_home, i) for i in negatives])
        else:
            true_negatives = self.get_fasta_ids(os.path.join(experiment_home, negatives))
            
        if type(positives) is list:
            true_positives = self.get_fasta_ids(
                [os.path.join(experiment_home, i) for i in positives])
        else:
            true_positives = self.get_fasta_ids(os.path.join(experiment_home, positives))
        
        target_psms = pd.read_csv(os.path.join(
            experiment_home, "percolator.target.psms.txt"), sep="\t")
        decoy_psms = pd.read_csv(os.path.join(
            experiment_home, "percolator.decoy.psms.txt"), sep="\t")

        return true_positives, true_negatives, target_psms, decoy_psms
