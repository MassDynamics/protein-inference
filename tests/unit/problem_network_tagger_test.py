import unittest
from protein_inference.inference.uniqueness_tagger import UniquenessTagger
from protein_inference.problem_network import ProblemNetwork
import networkx as nx


class ProblemNetworkTaggerTest(unittest.TestCase):

    def test_non_unique_tagging_true_and_false(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0)
        g.add_nodes_from([4,5], protein = 1)
        g.add_edges_from([(1,4),(2,4),(2,5),(3,5)])

        g = UniquenessTagger()._tag_non_unique_peptides(g)

        self.assertEqual(g.nodes[1]['unique'],1)
        self.assertEqual(g.nodes[2]['unique'],0)
        self.assertEqual(g.nodes[3]['unique'],1)

    def test_unique_evidenced_protein_true(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0)
        g.add_nodes_from([4,5], protein = 1)
        g.add_edges_from([(1,4),(2,4),(2,5),(3,5)])

        g = UniquenessTagger()._tag_non_unique_peptides(g)
        g = UniquenessTagger()._tag_unique_evidenced_protein(g)

        self.assertEqual(g.nodes[4]["unique_evidence"],1)
        self.assertEqual(g.nodes[5]["unique_evidence"],1)

    def test_unique_evidenced_protein_true_and_false(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0)
        g.add_nodes_from([4,5], protein = 1)
        g.add_edges_from([(1,4),(1,5),(2,4),(2,5),(3,5)])

        g = UniquenessTagger()._tag_non_unique_peptides(g)
        g = UniquenessTagger()._tag_unique_evidenced_protein(g)

        self.assertEqual(g.nodes[4]["unique_evidence"],0)
        self.assertEqual(g.nodes[5]["unique_evidence"],1)
