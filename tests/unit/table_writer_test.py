import unittest
from protein_inference.problem_network import ProblemNetwork
from protein_inference.table_maker import TableMaker
import networkx as nx
from copy import deepcopy


class TableMakerTest(unittest.TestCase):

    def test__get_edge_list_molecule_assignment(self):

        g = nx.Graph()
        g.add_nodes_from([1, 2, 3], protein=0)
        g.add_nodes_from([4, 5], protein=1)
        g.add_edges_from([(1, 4), (2, 4), (2, 5), (3, 5)], score=1)
        pn = ProblemNetwork(g)

        df = TableMaker()._get_edge_list(pn)

        for protein in pn.get_proteins():
            self.assertIn(protein, df["protein_id"].to_list())
            self.assertNotIn(protein, df["sequence_modified"].to_list())

        for peptide in pn.get_peptides():
            self.assertNotIn(peptide, df["protein_id"].to_list())
            self.assertIn(peptide, df["sequence_modified"].to_list())

    def test__get_edge_list_size(self):

        g = nx.Graph()
        g.add_nodes_from([1, 2, 3], protein=0)
        g.add_nodes_from([4, 5], protein=1)
        g.add_edges_from([(1, 4), (2, 4), (2, 5), (3, 5)], score=1)
        pn = ProblemNetwork(g)

        df = TableMaker()._get_edge_list(pn)

        self.assertEqual(df.shape[0], (4))
        self.assertEqual(df.shape[1], 3)

    def test_indistinguishable_col_basic(self):

        g = nx.Graph()
        g.add_nodes_from([1, 2, 3], protein=0)
        g.add_nodes_from([4, 5], protein=1, indistinguishable=[4, 5])
        g.add_edges_from([(1, 4), (2, 4), (2, 5), (3, 5)], score=1)
        pn = ProblemNetwork(g)

        df = TableMaker()._get_edge_list(pn)
        new_df = TableMaker().add_indistinguishable_col(pn, df)

        self.assertIn("indistinguishable", new_df.columns)

    def test_indistinguishable_col_assignment(self):

        g = nx.Graph()
        g.add_nodes_from([1, 2, 3], protein=0)
        g.add_nodes_from([4, 5], protein=1, indistinguishable=[4, 5])
        g.add_edges_from([(1, 4), (2, 4), (2, 5), (3, 5)], score=1)
        pn = ProblemNetwork(g)

        df = TableMaker()._get_edge_list(pn)
        new_df = TableMaker().add_indistinguishable_col(pn, df)

        self.assertEqual(new_df["indistinguishable"][0], [4, 5])
        self.assertEqual(new_df["indistinguishable"][1], [4, 5])

    def test_add_subset_col_basic(self):

        g = nx.Graph()
        g.add_nodes_from([1, 2, 3], protein=0)
        g.add_nodes_from([4, 5], protein=1, Group=4, major=4)
        g.add_edges_from([(1, 4), (2, 4), (2, 5), (3, 5)], score=1)
        pn = ProblemNetwork(g)

        df = TableMaker()._get_edge_list(pn)
        new_df = TableMaker().add_subset_col(pn, df)

        self.assertIn("subset", new_df.columns)

    def test_add_subset_col_assignment(self):
        '''not fully happy with this. Come back aftter editing add_subset col'''
        g = nx.Graph()
        g.add_nodes_from([1, 2, 3], protein=0)
        g.add_nodes_from([4, 5], protein=1, Group=4, major=4)
        g.add_edges_from([(1, 4), (2, 4), (2, 5), (3, 5)], score=1)
        pn = ProblemNetwork(g)

        df = TableMaker()._get_edge_list(pn)
        new_df = TableMaker().add_subset_col(pn, df)
        self.assertEqual(new_df["subset"][0], [4, 5])
        self.assertEqual(new_df["subset"][1], [4, 5])

    def test_get_protein_table_columns(self):

        g = nx.Graph()
        g.add_nodes_from([1, 3], protein=0, unique=1, razor=1)
        g.add_nodes_from([2], protein=0, unique=0)
        g.add_nodes_from([4], protein=1, score=1)
        g.add_nodes_from([5], protein=1, score=1)
        g.add_edges_from([(1, 4), (2, 4), (2, 5), (3, 5)], score=1)

        pn = ProblemNetwork(g)

        df = TableMaker().get_protein_table(pn)

        self.assertIn("score", df.columns)
        self.assertIn("unique", df.columns)
        self.assertIn("razor", df.columns)
        self.assertIn("non_unique", df.columns)
        self.assertIn("protein_id", df.columns)
        self.assertIn("sequence_modified", df.columns)

    def test_get_peptide_table_columns(self):

        g = nx.Graph()
        g.add_nodes_from([1], protein=0, unique=1,
                         razor=1, unique_evidence=True,
                         PEP = 0.01, q_value = 1)
        g.add_nodes_from([3], protein=0, unique=1,
                         razor=1, unique_evidence=True,
                         PEP = 0.4, q_value = 2)
        g.add_nodes_from([2], protein=0, unique=0, unique_evidence=False)
        g.add_nodes_from([4], protein=1, score=1, major=4)
        g.add_nodes_from([5], protein=1, score=1, major=4)
        g.add_edges_from([(1, 4), (2, 4), (2, 5), (3, 5)], score=1)

        pn = ProblemNetwork(g)

        df = TableMaker().get_peptide_table(pn)

        self.assertIn("score", df.columns)
        self.assertIn("unique", df.columns)
        self.assertIn("razor", df.columns)
        self.assertIn("major", df.columns)
        self.assertIn("protein_score", df.columns)
        self.assertIn("sequence", df.columns)
        self.assertIn("unique_evidence", df.columns)
