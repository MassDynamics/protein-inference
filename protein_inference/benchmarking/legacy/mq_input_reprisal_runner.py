from protein_inference.reprisal.greedy_algorithm import GreedyAlgorithm
from protein_inference.processing.psms_preprocessor import PSMsPreprocessor
from protein_inference.processing.psms_network_splitter import PSMsNetworkSplitter
from protein_inference.processing.psms_network_generator import PSMsNetworkGenerator
from protein_inference.inference.uniqueness_tagger import UniquenessTagger
from protein_inference.reprisal.greedy_algorithm import GreedyAlgorithm
from protein_inference.table_maker import TableMaker
from protein_inference.benchmarking.process_maxquant_output import ProcessMaxQuant
from protein_inference.inference.false_discovery_rate_calculator import FalseDiscoveryRateCalculator


class MQInputReprisalRunner():

    def run(self, mq_psms_table):

        target_psms = ProcessMaxQuant().make_psms_standard(mq_psms_table)
        #decoy_psms = ProcessMaxQuant().make_psms_standard(mq_psms_table, decoy=True)
        
        #might be able to merge all this with conventional runner, split over a couple difs
        target_psms = PSMsPreprocessor(target_psms, decoy=0).get_processed_psms()
        #decoy_psms = PSMsPreprocessor(decoy_psms, decoy=0).get_processed_psms()

        target_protein_table, target_peptide_table, solved_networks = self.get_output(target_psms)
        #decoy_protein_table, _, _ = self.get_output(decoy_psms)
#
        #target_protein_table = FalseDiscoveryRateCalculator().tag_FDR(target_protein_table,
        #                                                              decoy_protein_table)

        target_protein_table = target_protein_table.sort_values(
            ["score", "protein_id"], ascending=[False, True])

        #target_peptide_table = target_peptide_table.sort_values(
        #    ["score", "major", "sequence_modified"], ascending=[False, True, True])

        
        return target_protein_table, target_peptide_table, solved_networks

    def get_output(self,psms):  # a copy fron runner with merged removed

        network = PSMsNetworkGenerator(psms).generate_network()
        problem_networks = PSMsNetworkSplitter(network).split_networks()
        tagged_networks = UniquenessTagger().run(problem_networks)
        solved_networks = GreedyAlgorithm().run_system_wide(tagged_networks)
        #merged_networks = ProteinMerger().run_system_wide(solved_networks)
        protein_table = TableMaker().get_system_protein_table(solved_networks)
        peptide_table = TableMaker().get_system_peptide_table(solved_networks)


        return protein_table, peptide_table, solved_networks
