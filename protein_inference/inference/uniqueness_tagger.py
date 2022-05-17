from networkx import set_node_attributes
from protein_inference.problem_network import ProblemNetwork
from multiprocessing import cpu_count, Pool

class UniquenessTagger():

    def __init__(self):

        return

    def run(self, pn):

        network = self._tag_non_unique_peptides(pn.network)
        network = self._tag_unique_evidenced_protein(network)

        return ProblemNetwork(network)

    def _tag_non_unique_peptides(self, network):

        pn = ProblemNetwork(network)
        peptides = pn.get_peptides()
        list_uniqueness = list([1 == len(list(network.neighbors(peptide))) for
                                peptide in peptides])
        dict_uniqueness = dict(zip(peptides, list_uniqueness))

        # tag peptides with uniqueness
        set_node_attributes(network, dict_uniqueness, "unique")

        return network

    def _tag_unique_evidenced_protein(self, network):

        # empty dict of proteins inferred
        unique_evidence = {}

        pn = ProblemNetwork(network)

        for protein in pn.get_proteins():
            # calculates number of unique peptide neighbours
            n_unique_pep_neighbours = 0
            for n, d in network.nodes(data=True):
                if n in network.neighbors(protein):
                    if d['unique'] == 1:
                        n_unique_pep_neighbours = n_unique_pep_neighbours + 1

            if n_unique_pep_neighbours > 0:
                unique_evidence[protein] = True
            else:
                unique_evidence[protein] = False


        set_node_attributes(network, unique_evidence, "unique_evidence")

        return network



