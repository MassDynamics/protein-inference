import unittest
from protein_inference.problem_network import ProblemNetwork
from protein_inference.reprisal.greedy_algorithm import GreedyAlgorithm
import networkx as nx
from copy import deepcopy

class GreedyAlgorithmTest(unittest.TestCase):


    def test_score_protein_basic(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0, allocated = 0, unique = 0 )
        g.add_nodes_from([4,5], protein = 1)
        g.add_edges_from([(1,4),(2,4),(2,5),(3,5)], score = 1)
        pn = ProblemNetwork(g)

        score = GreedyAlgorithm().score_protein(pn, 4)

        self.assertEqual(score,2)

    def test_score_protein_allocated(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0, allocated = 1, unique = 0 )
        g.add_nodes_from([4,5], protein = 1)
        g.add_edges_from([(1,4),(2,4),(2,5),(3,5)], score = 1)
        pn = ProblemNetwork(g)

        score = GreedyAlgorithm().score_protein(pn, 4)

        self.assertEqual(score,0)

    def test_score_all_proteins_basic(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0, allocated = 0, unique = 0 )
        g.add_nodes_from([4,5], protein = 1)
        g.add_edges_from([(1,4),(2,4),(2,5),(3,5)], score = 1)
        pn = ProblemNetwork(g)

        scores = GreedyAlgorithm().score_all_proteins(pn)

        self.assertEqual(scores[4],2)
        self.assertEqual(scores[5],2)

    def test_score_all_proteins_asymmetric(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0, allocated = 0, unique = 0 )
        g.add_nodes_from([4,5], protein = 1)
        g.add_edges_from([(1,4),(2,5),(3,5)], score = 1)
        pn = ProblemNetwork(g)

        scores = GreedyAlgorithm().score_all_proteins(pn)

        self.assertEqual(scores[4],1)
        self.assertEqual(scores[5],2)

    def test_score_all_proteins_allocated(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0, allocated = 1, unique = 0 )
        g.add_nodes_from([4,5], protein = 1)
        g.add_edges_from([(1,4),(2,4),(2,5),(3,5)], score = 1)
        pn = ProblemNetwork(g)

        scores = GreedyAlgorithm().score_all_proteins(pn)

        self.assertEqual(scores[4],0)
        self.assertEqual(scores[5],0)


    def test_highest_scoring_proteins_basic(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0, allocated = 0, unique = 0 )
        g.add_nodes_from([4,5], protein = 1)
        g.add_edges_from([(1,4),(2,5),(3,5)], score = 1)
        pn = ProblemNetwork(g)

        protein = GreedyAlgorithm().get_highest_scoring_protein(pn)

        self.assertEqual(protein,5)

    def test_highest_scoring_proteins_tie(self):

        g = nx.Graph()
        g.add_nodes_from([1,2,3], protein = 0, allocated = 0, unique = 0 )
        g.add_nodes_from([4,5], protein = 1)
        g.add_edges_from([(1,4),(2,4),(2,5),(3,5)], score = 1)
        pn = ProblemNetwork(g)

        protein = GreedyAlgorithm().get_highest_scoring_protein(pn)

        self.assertEqual(protein,4)

    def test_highest_scoring_proteins_unique_only(self):

        g = nx.Graph()
        g.add_nodes_from([1], protein = 0, allocated = 0, unique = 1)
        g.add_nodes_from([2], protein = 0, allocated = 0, unique = 0)
        g.add_nodes_from([4], protein = 1, unique_evidence = True)
        g.add_nodes_from([5], protein = 1, unique_evidence = False)
        g.add_edges_from([(1,4),(2,4)], score = 1)
        g.add_edges_from([(2,5)], score = 10) #check if we still choose 4
        pn = ProblemNetwork(g)

        protein = GreedyAlgorithm().get_highest_scoring_protein(pn, unique_only=True)

        self.assertEqual(protein,4)

    def test_run_scores(self):

        g = nx.Graph()
        g.add_nodes_from([1], protein = 0, allocated = 0, unique = 1)
        g.add_nodes_from([2], protein = 0, allocated = 0, unique = 0)
        g.add_nodes_from([4], protein = 1, unique_evidence = True)
        g.add_nodes_from([5], protein = 1, unique_evidence = False)
        g.add_edges_from([(1,4),(2,4)], score = 1)
        g.add_edges_from([(2,5)], score = 10) #check if we still choose 4
        pn = ProblemNetwork(g)

        pn = GreedyAlgorithm().run(pn)

        self.assertEqual(pn.network.nodes[4]["score"], 2)
        self.assertEqual(pn.network.nodes[5]["score"], 0)

    def test_run_protein_allocator(self):

        g = nx.Graph()
        g.add_nodes_from([1], protein = 0, allocated = 0, unique = 1)
        g.add_nodes_from([2], protein = 0, allocated = 0, unique = 0)
        g.add_nodes_from([4], protein = 1, unique_evidence = True)
        g.add_nodes_from([5], protein = 1, unique_evidence = False)
        g.add_edges_from([(1,4),(2,4)], score = 1)
        g.add_edges_from([(2,5)], score = 10) #check if we still choose 4
        pn = ProblemNetwork(g)

        pn = GreedyAlgorithm().run(pn)

        self.assertEqual(pn.network.nodes[4]["major"], 4)
        self.assertEqual(pn.network.nodes[5]["major"], 4)

    def test_run_peptide_allocator(self):

        g = nx.Graph()
        g.add_nodes_from([1], protein = 0, allocated = 0, unique = 1)
        g.add_nodes_from([2], protein = 0, allocated = 0, unique = 0)
        g.add_nodes_from([4], protein = 1, unique_evidence = True)
        g.add_nodes_from([5], protein = 1, unique_evidence = False)
        g.add_edges_from([(1,4),(2,4)], score = 1)
        g.add_edges_from([(2,5)], score = 10) #check if we still choose 4
        pn = ProblemNetwork(g)

        pn = GreedyAlgorithm().run(pn)

        self.assertEqual(pn.network.nodes[1]["allocated"], 4)
        self.assertEqual(pn.network.nodes[2]["allocated"], 4)

    def test_run_razor_tagging(self):
        '''could possibly test this directly but it fits in well here'''

        g = nx.Graph()
        g.add_nodes_from([1], protein = 0, allocated = 0, unique = 1)
        g.add_nodes_from([2], protein = 0, allocated = 0, unique = 0)
        g.add_nodes_from([4], protein = 1, unique_evidence = True)
        g.add_nodes_from([5], protein = 1, unique_evidence = False)
        g.add_edges_from([(1,4),(2,4)], score = 1)
        g.add_edges_from([(2,5)], score = 10) #check if we still choose 4
        pn = ProblemNetwork(g)

        pn = GreedyAlgorithm().run(pn)

        self.assertEqual(pn.network.nodes[1]["razor"], True)
        self.assertEqual(pn.network.nodes[2]["razor"], True)

    def test_run_isomorphic(self):

        g = nx.Graph()
        g.add_nodes_from([1], protein = 0, allocated = 0, unique = 1)
        g.add_nodes_from([2], protein = 0, allocated = 0, unique = 0)
        g.add_nodes_from([4], protein = 1, unique_evidence = True)
        g.add_nodes_from([5], protein = 1, unique_evidence = False)
        g.add_edges_from([(1,4),(2,4)], score = 1)
        g.add_edges_from([(2,5)], score = 10) #check if we still choose 4
        
        pn = ProblemNetwork(g)
        pn2 = deepcopy(pn)

        pn = GreedyAlgorithm().run(pn)
        pn2 = GreedyAlgorithm().run(pn2)

        self.assertTrue(nx.is_isomorphic(pn.network,pn2.network))

    def test_run_systemwide(self):

        g = nx.Graph()
        g.add_nodes_from([1], protein = 0, allocated = 0, unique = 1)
        g.add_nodes_from([2], protein = 0, allocated = 0, unique = 0)
        g.add_nodes_from([4], protein = 1, unique_evidence = True)
        g.add_nodes_from([5], protein = 1, unique_evidence = False)
        g.add_edges_from([(1,4),(2,4)], score = 1)
        g.add_edges_from([(2,5)], score = 10) #check if we still choose 4
        
        pn = ProblemNetwork(g)
        pn2 = deepcopy(pn)

        pns = GreedyAlgorithm().run_system_wide([pn,pn2])

        pn1_non_par = GreedyAlgorithm().run(pn)
        pn2_non_par = GreedyAlgorithm().run(pn2)
        self.assertTrue(nx.is_isomorphic(pns[0].network, pn1_non_par.network))
        self.assertTrue(nx.is_isomorphic(pns[1].network, pn2_non_par.network))