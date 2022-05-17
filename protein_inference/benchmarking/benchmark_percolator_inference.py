import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


class ProcessPercolatorInference():

    def load_protein_table(self, path):
        protein_table = pd.read_csv(path, sep="\t")
        protein_table["peptides"] = protein_table.peptideIds.str.split(r"\s")
        return protein_table[["ProteinId", "ProteinGroupId", "q-value","posterior_error_prob", "peptides"]]

    def label_problem_networks(self, pns, percolator_protein_table):
        '''Since the original workflow is built around percolator output,
        we can assume the pn's are accurate and find proteins in the generated
        network and assign the scores and groups according to this table'''

    def label_problem_network(self, network, percolator_protein_table):
        '''Since the original workflow is built around percolator output,
        we can assume the pn's are accurate and find proteins in the generated
        network and assign the scores and groups according to this table'''

        for _, row in percolator_protein_table.iterrows():
            # get values
            protein = row["ProteinId"]
            group = row["ProteinGroupId"]
            q_value = row["q-value"]
            PEP = row["posterior_error_prob"]
            peptides = row["peptides"]

            # assign
            if protein in network.nodes():
                network.nodes[protein]["major"] = group
                network.nodes[protein]["q_value"] = q_value
                network.nodes[protein]["PEP"] = PEP

            for peptide in peptides:
                if peptide in network.nodes():
                    network.nodes[peptide]["allocated"] = group
                # There are probably several good reasons peptides are lost
                # in our network table during preprocessing.
                # except:
                #    print("peptide", peptide, " is not in the table")

        return network


class EvaluatePercolatorInference():
    '''This is only about comparign q-value's and FDR scores.
    Code has been written with expectation of ground truth data in mind.'''

    def plot_percolator_reprisal_predictions(self, reprisal_PI_df, percolator_PI_df,
                                             real=False, fake=False):
        comparison_df = pd.DataFrame.merge(reprisal_PI_df, percolator_PI_df,
                                           how="outer",
                                           left_on="protein_id",
                                           right_on="ProteinId",
                                           suffixes=["_rep", "_per"])

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines",
                                 line=go.scatter.Line(
                                     color="gray", dash="dashdot"),
                                 showlegend=False))
        fig.add_trace(go.Scatter(x=[0, 1], y=[0.01, 0.01], mode="lines",
                                 line=go.scatter.Line(
                                     color="gray", dash="dashdot"),
                                 showlegend=False))
        fig.add_trace(go.Scatter(x=[0.01, 0.01], y=[0, 1], mode="lines",
                                 line=go.scatter.Line(
                                     color="gray", dash="dashdot"),
                                 showlegend=False))

        if real and fake:
            tmp = comparison_df[comparison_df.protein_id.apply(
                lambda x: x in fake)]

            fig.add_trace(go.Scatter(x=tmp["q-value_rep"], y=tmp["q-value_per"], mode='markers',
                                     line=dict(color="DarkRed"), 
                                     name="Entrapment (Doesn't exist)", hovertext=tmp.ProteinId))
            tmp = comparison_df[comparison_df.protein_id.apply(
                lambda x: x in real)]
            fig.add_trace(go.Scatter(x=tmp["q-value_rep"], y=tmp["q-value_per"], mode='markers',
                                     line=dict(color="LightSeaGreen"), 
                                     name="True Target (exists)", hovertext=tmp.ProteinId))
        else:
            tmp = comparison_df
            fig.add_trace(go.Scatter(x=tmp["q-value_rep"], y=tmp["q-value_per"], mode='markers',
                                     line=dict(color="DarkRed"), name="All (no ground truth give)", 
                                     hovertext=tmp.ProteinId))
        
        fig.update_layout(title="REPRISAL - PERCOLATOR Benchmark",
                            xaxis_title="REPRISAL q-value",
                            yaxis_title="Percolator q-value")

        
        fig.update_xaxes(range=[0, 0.2])
        fig.update_yaxes(range=[0, 0.2])
        fig.show()
        return

    def compare_perc_rep_fdr(self, perc_protein_table, rep_protein_table):
        comparison_df = pd.DataFrame.merge(rep_protein_table, perc_protein_table,
                                           how="outer",
                                           left_on="protein_id",
                                           right_on="ProteinId",
                                           suffixes=["_rep", "_per"])

        tmp = comparison_df.rename(columns={"FDR": "REPRISAL FDR",
                                            "q-value": "PERCOLATOR q-value"})
        fig = px.scatter(tmp, x="REPRISAL FDR", y="PERCOLATOR q-value",
                         hover_data=["ProteinId", "REPRISAL FDR", "PERCOLATOR q-value"])
        fig.show()

    def check_perc_rep_agreement(self, perc_protein_table, rep_protein_table):
        comparison_df = pd.DataFrame.merge(rep_protein_table, perc_protein_table,
                                           how="outer",
                                           left_on="protein_id",
                                           right_on="ProteinId",
                                           suffixes=["_rep", "_per"])

        rep_not_perc = comparison_df[(comparison_df.FDR < 0.01) & (
            comparison_df["q-value"] > 0.01)]
        perc_not_rep = comparison_df[(comparison_df.FDR > 0.01) & (
            comparison_df["q-value"] < 0.01)]

        return rep_not_perc, perc_not_rep
