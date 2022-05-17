import unittest
from protein_inference.problem_network import ProblemNetwork
from protein_inference.inference.protein_merger import ProteinMerger
import networkx as nx
from copy import deepcopy

class ProteinMergerTest(unittest.TestCase):

    def test_get_mapping_basic(self):

        g = nx.Graph()
        g.add_nodes_from(["1"], protein = 0)
        g.add_nodes_from(["2"], protein = 0)
        g.add_nodes_from(["4"], protein = 1)
        g.add_nodes_from(["5"], protein = 1)
        g.add_edges_from([("1","4"),("2","4")], score = 1)
        g.add_edges_from([("2","5")], score = 10) 
        pn = ProblemNetwork(g)

        df = ProteinMerger().get_mapping(pn)
        df = df.sort_values("protein")

        self.assertEqual(df["peptides"][0],["1","2"])
        self.assertEqual(df["peptides"][1],["2"])

    def test_get_mapping_indistinguishable(self):

        g = nx.Graph()
        g.add_nodes_from(["1"], protein = 0)
        g.add_nodes_from(["2"], protein = 0)
        g.add_nodes_from(["4"], protein = 1)
        g.add_nodes_from(["5"], protein = 1)
        g.add_edges_from([("1","4"),("2","4")], score = 1)
        g.add_edges_from([("2","5"),("1","5")], score = 10) 
        pn = ProblemNetwork(g)

        df = ProteinMerger().get_mapping(pn)
        df = df.sort_values("protein")

        self.assertEqual(df["peptides"][0],["1","2"])
        self.assertEqual(df["protein"][0],["4","5"])

    def test_get_named_proteins_basic(self):

        g = nx.Graph()
        g.add_nodes_from(["1"], protein = 0)
        g.add_nodes_from(["2"], protein = 0)
        g.add_nodes_from(["4"], protein = 1)
        g.add_nodes_from(["5"], protein = 1)
        g.add_edges_from([("1","4"),("2","4")], score = 1)
        g.add_edges_from([("2","5")], score = 10) 
        g.nodes["4"]["score"] = 2
        g.nodes["5"]["score"] = 0

        pn = ProblemNetwork(g)

        df = ProteinMerger().get_named_proteins(pn)
        df = df.sort_values("protein")

        self.assertEqual(df["named"][0],"4")
        self.assertEqual(df["named"][1],"5")

    def test_get_named_proteins_indistinguishable(self):

        g = nx.Graph()
        g.add_nodes_from(["1"], protein = 0)
        g.add_nodes_from(["2"], protein = 0)
        g.add_nodes_from(["4"], protein = 1)
        g.add_nodes_from(["5"], protein = 1)
        g.add_edges_from([("1","4"),("2","4")], score = 1)
        g.add_edges_from([("2","5"),("1","5")], score = 10) 
        g.nodes["4"]["score"] = 2
        g.nodes["5"]["score"] = 0

        pn = ProblemNetwork(g)

        df = ProteinMerger().get_named_proteins(pn)
        df = df.sort_values("protein")

        self.assertEqual(df["named"][0],"4")

    def test_get_named_proteins_indistinguishable_tie(self):

        g = nx.Graph()
        g.add_nodes_from(["1"], protein = 0)
        g.add_nodes_from(["2"], protein = 0)
        g.add_nodes_from(["4"], protein = 1)
        g.add_nodes_from(["5"], protein = 1)
        g.add_edges_from([("1","4"),("2","4")], score = 1)
        g.add_edges_from([("2","5"),("1","5")], score = 10) 
        g.nodes["4"]["score"] = 2
        g.nodes["5"]["score"] = 2

        pn = ProblemNetwork(g)

        df = ProteinMerger().get_named_proteins(pn)
        df = df.sort_values("protein")

        self.assertEqual(df["named"][0],"4")

    def test_run_no_merges(self):

        g = nx.Graph()
        g.add_nodes_from(["1"], protein = 0)
        g.add_nodes_from(["2"], protein = 0)
        g.add_nodes_from(["4"], protein = 1)
        g.add_nodes_from(["5"], protein = 1)
        g.add_edges_from([("1","4"),("2","4")], score = 1)
        g.add_edges_from([("2","5")], score = 10) 
        g.nodes["4"]["score"] = 2
        g.nodes["5"]["score"] = 0

        pn = ProblemNetwork(g)
        pn = ProteinMerger().run(pn)
        
        self.assertEqual(len(pn.get_proteins()),2)

    def test_run_a_merge(self):

        g = nx.Graph()
        g.add_nodes_from(["1"], protein = 0)
        g.add_nodes_from(["2"], protein = 0)
        g.add_nodes_from(["4"], protein = 1)
        g.add_nodes_from(["5"], protein = 1)
        g.add_edges_from([("1","4"),("2","4")], score = 1)
        g.add_edges_from([("2","5"),("1","5")], score = 10) 
        g.nodes["4"]["score"] = 2
        g.nodes["5"]["score"] = 0

        pn = ProblemNetwork(g)
        pn = ProteinMerger().run(pn)
        
        self.assertEqual(len(pn.get_proteins()),1)
        self.assertEqual(pn.get_proteins(),["4"])

    def test_run_indistinguishable_label(self):

        g = nx.Graph()
        g.add_nodes_from(["1"], protein = 0)
        g.add_nodes_from(["2"], protein = 0)
        g.add_nodes_from(["4"], protein = 1)
        g.add_nodes_from(["5"], protein = 1)
        g.add_edges_from([("1","4"),("2","4")], score = 1)
        g.add_edges_from([("2","5"),("1","5")], score = 10) 
        g.nodes["4"]["score"] = 2
        g.nodes["5"]["score"] = 0

        pn = ProblemNetwork(g)
        pn = ProteinMerger().run(pn)
        
        self.assertEqual(pn.network.nodes["4"]["indistinguishable"],["5"])

    def test_run_isomorphic(self):
        
        g = nx.Graph()
        g.add_nodes_from(["1"], protein = 0)
        g.add_nodes_from(["2"], protein = 0)
        g.add_nodes_from(["4"], protein = 1)
        g.add_nodes_from(["5"], protein = 1)
        g.add_edges_from([("1","4"),("2","4")], score = 1)
        g.add_edges_from([("2","5"),("1","5")], score = 10) 
        g.nodes["4"]["score"] = 2
        g.nodes["5"]["score"] = 0

        pn = ProblemNetwork(g)
        pn2 = deepcopy(pn)


        pn1 =  ProteinMerger().run(pn)
        pn2 =  ProteinMerger().run(pn2)

        self.assertTrue(nx.is_isomorphic(pn1.network,pn2.network))

    def test_run_system_wide(self):

        g = nx.Graph()
        g.add_nodes_from(["1"], protein = 0)
        g.add_nodes_from(["2"], protein = 0)
        g.add_nodes_from(["4"], protein = 1)
        g.add_nodes_from(["5"], protein = 1)
        g.add_edges_from([("1","4"),("2","4")], score = 1)
        g.add_edges_from([("2","5"),("1","5")], score = 10) 
        g.nodes["4"]["score"] = 2
        g.nodes["5"]["score"] = 0

        pn = ProblemNetwork(g)
        pn2 = deepcopy(pn)

        pns = []
        pns.append(ProteinMerger().run(pn))
        pns.append(ProteinMerger().run(pn2))

        pn1_non_par =  ProteinMerger().run(pn)
        pn2_non_par =  ProteinMerger().run(pn2)
        self.assertTrue(nx.is_isomorphic(pns[0].network, pn1_non_par.network))
        self.assertTrue(nx.is_isomorphic(pns[1].network, pn2_non_par.network))