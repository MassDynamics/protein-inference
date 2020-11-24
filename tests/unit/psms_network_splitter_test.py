import unittest
from protein_inference.problem_network import ProblemNetwork
from protein_inference.processing.psms_network_splitter import PSMsNetworkSplitter
import networkx as nx

class PSMsNetworkSplitterTest(unittest.TestCase):

    def test_split_networks_no_splitting(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0)
        g.add_nodes_from([4,5], protein = 1)
        g.add_edges_from([(1,4),(2,4),(2,5),(3,5)])

        list_subnetworks = PSMsNetworkSplitter(g).split_networks()

        self.assertEqual(len(list_subnetworks), 1)

    def test_split_networks_some_splitting(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0)
        g.add_nodes_from([4,5], protein = 1)
        g.add_edges_from([(1,4),(2,4),(3,5)])

        list_subnetworks = PSMsNetworkSplitter(g).split_networks()

        self.assertEqual(len(list_subnetworks), 2)