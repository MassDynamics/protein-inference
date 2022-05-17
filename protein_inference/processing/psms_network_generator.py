from networkx import Graph

class PSMsNetworkGenerator:
    """
    A class that generates a networkx graph from a psm
    table.

    Attributes:
    -----------
    PSMs : ProcessedPSMs
        A ProcessedPSMs instance. 

    Methods:
    --------
    generate_network()
        Generates the networkx graph with appropriate node
        types, labels and edge attributes.
    
    """

    def __init__(self, PSMs):

        self.PSMs = PSMs

        return

    def generate_network(self):

        network = Graph()

        network.add_nodes_from(self.PSMs.get_peptides(),
                               protein=0)
        network.add_nodes_from(self.PSMs.get_proteins(),
                               protein=1)
        # add edges
        for _, row in self.PSMs.df.iterrows():
            network.add_edge(row["protein id"],
                             row["sequence"],
                             PEP=row["PEP"],
                             score=row["score"],
                             q_value=row["q_value"])

        return network
