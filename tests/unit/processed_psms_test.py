import unittest
from protein_inference.processing.processed_psms import ProcessedPSMs
import pandas as pd

class ProcessedPSMsTest(unittest.TestCase):

    def test_init_df_type(self):

        df = pd.read_csv("tests/test_data/test_psm.txt", sep = "\t")
        pro_psms = ProcessedPSMs(df)
        
        self.assertEqual(type(pro_psms.df), type(pd.DataFrame()))


    def test_get_proteins_complete(self):

        df = pd.read_csv("tests/test_data/test_psm.txt", sep = "\t")
        pro_psms = ProcessedPSMs(df)
        
        self.assertIn("sp|P25911|LYN_MOUSE", pro_psms.get_proteins())
        self.assertIn("sp|P47955|RLA1_MOUSE", pro_psms.get_proteins())

    def test_get_proteins_size(self):

        df = pd.read_csv("tests/test_data/test_psm.txt", sep = "\t")
        pro_psms = ProcessedPSMs(df)
        
        self.assertEqual(2, len(pro_psms.get_proteins()))


    def test_get_peptides_complete(self):
        df = pd.read_csv("tests/test_data/test_psm.txt", sep = "\t")
        pro_psms = ProcessedPSMs(df)
        
        self.assertIn("ALANVNIGSLICNVGAGGPAPAAGAAPAGGAAPSTAAAPAEEK", 
                        pro_psms.get_peptides())
        self.assertIn("DPEEQGDIVVALYPYDGIHPDDLSFK", 
                        pro_psms.get_peptides())

    def test_get_peptides_size(self):

        df = pd.read_csv("tests/test_data/test_psm.txt", sep = "\t")
        pro_psms = ProcessedPSMs(df)
        
        self.assertEqual(2, len(pro_psms.get_peptides()))

    def test_get_shape(self):
        df = pd.read_csv("tests/test_data/test_psm.txt", sep = "\t")
        pro_psms = ProcessedPSMs(df)
        
        self.assertEqual((4,13),pro_psms.get_shape())


