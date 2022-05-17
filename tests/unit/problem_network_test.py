import unittest
from protein_inference.problem_network import ProblemNetwork
import networkx as nx

class ProblemNetworkTest(unittest.TestCase):

    def test_get_proteins_several(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 1)
        
        pn = ProblemNetwork(g)

        self.assertEqual(pn.get_proteins(),[1,2,3])

    def test_get_proteins_none(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0)
        
        pn = ProblemNetwork(g)

        self.assertEqual(pn.get_proteins(),[])


    def test_get_peptides_several(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0)
        g.add_nodes_from([4,5], protein = 1)

        pn = ProblemNetwork(g)

        self.assertEqual(pn.get_peptides(),[1,2,3])

    def test_update_nodes_no_prior_att(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 1)
        pn = ProblemNetwork(g)
        pn.update_nodes([1,2,3], "att",3)

        self.assertEqual(pn.network.nodes[1]["att"],3)

    def test_update_nodes_prior_atts(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 1)
        pn = ProblemNetwork(g)
        pn.update_nodes([1,2,3], "att",3)

        pn.update_nodes([2,3], "att",4)

        self.assertEqual(pn.network.nodes[1]["att"],3)
        self.assertEqual(pn.network.nodes[2]["att"],4)

    def test_pick_nodes_all(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 1)
        pn = ProblemNetwork(g)

        self.assertEqual(pn.pick_nodes("protein",1), [1,2,3])

    def test_pick_nodes_none(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 1)
        pn = ProblemNetwork(g)

        self.assertEqual(pn.pick_nodes("protein",0), [])

    def test_get_node_attribute_dict_all(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 1)
        pn = ProblemNetwork(g)

        self.assertEqual(pn.get_node_attribute_dict("protein"), {1:1,2:1,3:1})

    
    def test_get_node_attribute_dict_some(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 1)
        g.add_nodes_from([4,5])
        pn = ProblemNetwork(g)

        self.assertEqual(pn.get_node_attribute_dict("protein"), {1:1,2:1,3:1})

    def test_get_node_attribute_dict_none(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 1)
        pn = ProblemNetwork(g)

        self.assertEqual(pn.get_node_attribute_dict("peptide"), {})