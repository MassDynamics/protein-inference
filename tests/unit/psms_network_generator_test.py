import unittest
from protein_inference.processing.processed_psms import ProcessedPSMs
from protein_inference.processing.psms_network_generator import PSMsNetworkGenerator
import pandas as pd

class PSMsNetworkGeneratorTest(unittest.TestCase):

    def test_generate_network_nodes(self):

        #need psms object
        df = pd.read_csv("tests/test_data/test_psm.txt", sep = "\t")
        df = df.head(2)
        pro_psms = ProcessedPSMs(df)

        network = PSMsNetworkGenerator(pro_psms).generate_network()

    
        self.assertEqual(len(network.nodes), 4)


    def test_generate_network_edges(self):

        #need psms object
        df = pd.read_csv("tests/test_data/test_psm.txt", sep = "\t")
        pro_psms = ProcessedPSMs(df)

        network = PSMsNetworkGenerator(pro_psms).generate_network()

    
        self.assertEqual(len(network.edges), 2)

    def test_generate_network_edge_properties(self):

        #need psms object
        df = pd.read_csv("tests/test_data/test_psm.txt", sep = "\t")
        pro_psms = ProcessedPSMs(df)

        network = PSMsNetworkGenerator(pro_psms).generate_network()

        edge_data = network.get_edge_data("sp|P47955|RLA1_MOUSE",
        "ALANVNIGSLICNVGAGGPAPAAGAAPAGGAAPSTAAAPAEEK")

        self.assertIn("PEP", edge_data)
        self.assertIn("score", edge_data)
