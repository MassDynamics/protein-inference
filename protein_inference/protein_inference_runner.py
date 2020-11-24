import os
from os.path import join
from pandas import read_csv
import pickle
from multiprocessing import cpu_count, Pool

from protein_inference.processing.psms_preprocessor import PSMsPreprocessor
from protein_inference.processing.psms_network_generator import PSMsNetworkGenerator
from protein_inference.processing.psms_network_splitter import PSMsNetworkSplitter

from protein_inference.table_maker import TableMaker

from protein_inference.inference.uniqueness_tagger import UniquenessTagger
from protein_inference.inference.protein_merger import ProteinMerger
from protein_inference.inference.false_discovery_rate_calculator import FalseDiscoveryRateCalculator

from protein_inference.reprisal.greedy_algorithm import GreedyAlgorithm


class ProteinInferenceRunner():
    """
    The main class that drives the protein inference workflow. 

    This class groups the main coordinating functions that organise
    a protein inference workflow including PSM processing, graph
    generation, splitting, tagging, scoring, fdr calculations
    and output file writing. 

    To instantiate:

    >>> runner = ProteinInferenceRunner()
    """

    def run(self, target_path, decoy_path, output_directory, scoring_method=GreedyAlgorithm):
        """
        This method calls the protein inference workflow.

        This method locates two psm tables, corresponding to target and decoy matches,
        and uses them to generate output peptide and protein tables with scores and
        false discovery rates indicating which proteins are inferred and which
        peptides have been assigned to each protein (as evidence). 

        Parameters:
        -----------
        target_path : string
            A path to a valid psm table with target matches
        decoy_path : string
            A path to a valid psm table with decoy matches
        output_directory : string
            A path to the output directory
        scoring_method : scoring method object
            A valid scoring method such as those in inference.scorers

        Returns
        -------
        reprisal.target.proteins.csv : csv
            A file in the output directory. Describes the inferences made. 
        reprisal.target.peptides.csv : csv
            A file in the output directory. Describes the inferences made.    
        target_networks.p : pickle file
            A file that can be loaded using pickle to retrieve the target networks
            for visualization. 

        """

        print("Loading data...")
        target_psms = read_csv(target_path, sep="\t")
        decoy_psms = read_csv(decoy_path, sep="\t")

        print("Scoring Decoys...")
        decoy_output = self.get_output(
            decoy_psms, decoy=1, scoring_method=scoring_method)

        print("Scoring Targets...")
        target_output = self.get_output(
            target_psms, scoring_method=scoring_method)

        del(target_psms)
        del(decoy_psms)
        
        decoy_protein_table = decoy_output[0]
        decoy_peptide_table = decoy_output[1]
        del(decoy_output)

        target_protein_table = target_output[0]
        target_peptide_table = target_output[1]
        target_networks = target_output[2]
        del(target_output)
        
        print("Estimating False Discovery Rates...")
        target_fdr_table = FalseDiscoveryRateCalculator().tag_q_value(target_protein_table[target_protein_table.score > 0],
                                                                      decoy_protein_table[decoy_protein_table.score > 0])

        target_fdr_table = target_fdr_table.sort_values(
            ["q-value", "ProteinId"], ascending=[True, True])

        target_peptide_table = target_peptide_table.sort_values(
            ["score", "major", "sequence"], ascending=[False, True, True])
        
        print("Writing Outputs...")
        self.write_tables(target_fdr_table,
                          target_peptide_table, 
                          decoy_protein_table,
                          decoy_peptide_table,
                          output_directory)
        pickle.dump(target_networks, open(os.path.join(
            output_directory, "target_networks.p"), "wb"))

    def get_output(self, psms, decoy=0, scoring_method=GreedyAlgorithm):
        """
        Processes either target or decoy data. 
        
        Preprocesses data,generates, splits and tags networks. 
        Scores proteins, merges proteins with identical neighbors 
        and retrieves protein and peptide tables.

        Parameters
        ----------
        psms : pandas.DataFrame
            A pandas dataframe contain a psms table. 
        decoy: bool
            A boolean indicating whether decoys should be removed
            from the input table.
        scoring_method : scoring method object
            A valid scoring method such as those in inference.scorers

        Returns
        -------
        protein_table : pandas.DataFrame
            A pandas dataframe containing proteins and inference metrics.
        peptide_table : pandas.DataFrame
            A pandas dataframe containing peptides and assignments. 
        solved_networks: a list of problem networks with scores and node
            labels. 

        """
        psms = PSMsPreprocessor(psms, decoy=decoy).get_processed_psms()
        network = PSMsNetworkGenerator(psms).generate_network()
        del(psms)

        problem_networks = PSMsNetworkSplitter(network).split_networks()
        del(network)
        tagged_networks = self.parallel_apply(
            problem_networks, UniquenessTagger().run)
        del(problem_networks)

        solved_networks = self.parallel_apply(
            tagged_networks, scoring_method().run)
        del(tagged_networks)
        
        merged_networks = self.parallel_apply(
            solved_networks, ProteinMerger().run)
        
        protein_table = TableMaker().get_system_protein_table(merged_networks)
        peptide_table = TableMaker().get_system_peptide_table(merged_networks)

        return protein_table, peptide_table, solved_networks

    def write_tables(self, target_protein_table, target_peptide_table, decoy_protein_table, decoy_peptide_table, output_directory):
        """
        Writes tables to csvs.
        
        Parameters:
        -----------
        target_protein_table : pandas.DataFrame
            A pandas dataframe containing proteins and inference metrics.
        target_peptide_table : pandas.DataFrame
            A pandas dataframe containing peptides and assignments. 
        decoy_protein_table : pandas.DataFrame
            A pandas dataframe containing proteins and inference metrics.
        decoy_peptide_table : pandas.DataFrame
            A pandas dataframe containing peptides and assignments. 
        output_directory : string
            A path to the output directory

        Returns
        -------
        reprisal.target.proteins.csv : csv
            A file in the output directory. Describes the inferences made. 
        reprisal.target.peptides.csv : csv
            A file in the output directory. Describes the inferences made. 
        reprisal.decoy.proteins.csv : csv
            A file in the output directory. Describes the inferences made. 
        reprisal.decoy.peptides.csv : csv
            A file in the output directory. Describes the inferences made. 
        
        """

        target_protein_table.to_csv(
            join(output_directory, "reprisal.target.proteins.csv"), index=False)
        target_peptide_table.to_csv(
            join(output_directory, "reprisal.target.peptides.csv"), index=False)
        decoy_protein_table.to_csv(
            join(output_directory, "reprisal.decoy.proteins.csv"), index=False)
        decoy_peptide_table.to_csv(
            join(output_directory, "reprisal.decoy.peptides.csv"), index=False)
        return

    def parallel_apply(self, pns, func):
        '''
        Uses multiprocessing module in mython to apply a function
        to every problem network in a list of networks
        
        Parameters:
        -----------
        pns : list
            A list of ProblemNetworks.
        func: function or method
            A function or method that acts on a ProblemNetwork.
            
        Returns:
        --------
        pns : list
            A list of processed ProblemNetworks.
        '''
        p = Pool(cpu_count())
        pns = p.map(func, pns)

        return pns
