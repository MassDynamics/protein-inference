from networkx import set_node_attributes
import pprint
from itertools import combinations
from copy import deepcopy


class ProblemNetwork:
    '''

    Stores a networkx object with nodes for proteins and 
    peptides and edges for PSMs (protein-spectral-matches).


    This object stores the protein-peptide PSM network along
    with several handy methods that support tagging, scoring
    and plotting functionalities.

    Attributes
    ----------
    network : nx.Graph
        A networkx graph storing the ProblemNetwork

    '''

    def __init__(self, network):

        self.network = network

    def get_proteins(self):
        """Return a list of the proteins contained in the network."""
        return list([n for n, d in self.network.nodes(data=True) if d['protein'] == 1])

    def get_peptides(self):
        """Returns a list of peptides contained in the network."""
        return list([n for n, d in self.network.nodes(data=True) if d['protein'] == 0])

    def update_nodes(self, nodes, attribute, value):
        """
        
        Takes a list of nodes in the networks and gives them
        a node *attribute* with a particular *value*.
        
        Parameters
        ----------
        nodes : list
            A list of nodes ids (strings) to update.
        attribute: string
            The name of the attribute to update
        value : NA
            The value to set the attribute to. Can be string or float.
        """

        att_dict = dict(zip(nodes, [value]*len(nodes)))

        # update nodes
        set_node_attributes(self.network,
                            att_dict,
                            attribute)

        return

    def pick_nodes(self, attribute, value):
        """
        Returns a list of all nodes with a particular 
        attribute - value pairing
        
        Parameters
        ----------
        attribute: string
            The node attribute.
        value : NA
            The value the attribute must be set to.
        
        """


        picked = []
        for n, d in self.network.nodes(data=True):
            if attribute in d.keys():
                if d[attribute] == value:
                    picked.append(n)

        return picked

    def print_nodes(self):
        """Uses pretty print to show all nodes and their
        attributes."""
        pp = pprint.PrettyPrinter(indent=4,)
        pp.pprint(dict(self.network.nodes(data=True)))

        return

    def print_edges(self):
        """Uses pretty print to show all edges and their
        attributes."""
        pp = pprint.PrettyPrinter(indent=4,)
        edge_rep = {(u[0], u[1]): u[2] for u in self.network.edges(data=True)}
        pp.pprint(edge_rep)

    def get_node_attribute_dict(self, attribute):
        """
        Returns a dictionary of all node:value
        pairings for a particular value.
        
        Parameters
        ----------
        attribute: string
            The name of the attribute for which to retrieve the
            node:value pairings.
        """
        out = {}

        for node in self.network.nodes():
            if attribute in self.network.nodes[node].keys():
                out[node] = self.network.nodes[node][attribute]

        return out

    def get_edge_attribute_dict(self, attribute):
        """
        Returns a dictionary of all edge:value
        pairings for a particular value.
        
        Parameters
        ----------
        attribute: string
            The name of the attribute for which to retrieve the
            edge:value pairings.
        """
        out = {}

        for edge in self.network.edges():
            if attribute in self.network.edges[edge].keys():
                out[edge] = self.network.edges[edge][attribute]

        return out
