from .process_maxquant_output import ProcessMaxQuantOutput
from .process_protein_inference_output import ProcessProteinInferenceOutput
import matplotlib.pyplot as plt

class CompareMaxQuantOutput():
    '''This class collects functiosn for the comparison 
    of a maxquant protein groups file (read in via pandas)
    to the protein table output.'''

    def compare_infered_proteins(self, mqoutput_df, protein_table):
        '''This function summarizes the number of proteins inferred
        in the protein table and maxquant output file.'''

        out = {}
        set_mq = ProcessMaxQuantOutput().get_major_proteins(mqoutput_df)
        set_pi = ProcessProteinInferenceOutput().get_major_proteins(protein_table)
        out["n_pimd"] = len(set_pi)
        out["n_maxquant"] = len(set_mq)
        out["n_both"] = len(set_pi.intersection(set_mq))
        out["n_pimd_not_maxquant"] = len(set_pi.union(set_mq).difference(set_mq))
        out["n_maxquant_not_pimd"] = len(set_pi.union(set_mq).difference(set_pi))

        return out

    def compare_protein_grouping(self, mqoutput_df, protein_table):

        out = {}
        set_mq = ProcessMaxQuantOutput().get_protein_groups_dictionary(mqoutput_df)
        set_pi = ProcessProteinInferenceOutput().get_protein_groups_dictionary(protein_table)


        out["ave_group_size_mq"] = self.average_set_size(set_mq)
        out["ave_group_size_pimd"] = self.average_set_size(set_pi)
        #jaccard sim
        out["ave_jac_sim_shared_groups"] = self.get_average_jaccard_sim(set_mq, set_pi)

        return out


    def jaccard_sim(self,set_a,set_b):
        return len(set_a.intersection(set_b))/len(set_a.union(set_b))

    def get_average_jaccard_sim(self,dict_of_sets_a, dict_of_sets_b):

        common_groups = list([key for key in dict_of_sets_a.keys() if key in dict_of_sets_b.keys()])
        similarities = []
        for group in common_groups:
            similarities.append(self.jaccard_sim(dict_of_sets_a[group], 
                                        dict_of_sets_b[group]))
        return sum(similarities)/len(similarities)
    
    def average_set_size(self,dict_of_sets):

        lengths = [len(protein_set) for protein_set in dict_of_sets.values()]
        return sum(lengths)/len(lengths)
        
    def compare_protein_scores(self, mqoutput_df, protein_table):

        #in case not uniprot
        protein_table = ProcessProteinInferenceOutput().convert_protein_table_to_uniprot(protein_table)
        
        score_dict_mq = dict(zip(mqoutput_df["Majority protein IDs"], mqoutput_df.Score))
        score_dict_pimd = dict(zip(protein_table.protein_id, protein_table.score))

        common_groups = list([key for key in score_dict_mq.keys() if key in score_dict_pimd.keys()])
        
        protein_ids = []
        pimd_scores = []
        mq_scores = []

        for id in common_groups:
            protein_ids.append(id)
            pimd_scores.append(score_dict_pimd[id])
            mq_scores.append(score_dict_mq[id])

        plt.figure(figsize=[10,10])
        plt.scatter(pimd_scores, mq_scores)
        plt.xlabel("MD Protein Inference Scores")
        plt.ylabel("Max Quant Protein Inference Scores")
        plt.title("Score Comparison")
        plt.savefig("ScoreComparison.png")