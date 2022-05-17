from pandas import read_csv
from protein_inference.processing.psms_preprocessor import PSMsPreprocessor
from protein_inference.processing.psms_network_splitter import PSMsNetworkSplitter
from protein_inference.processing.psms_network_generator import PSMsNetworkGenerator

class ProcessMaxQuantOutput:

    def run(self,path):

        mqoutput_df = self.load(path)
        mqoutput_df = self._remove_decoys(mqoutput_df)

        return mqoutput_df

    def load(self, path):

        mqoutput_df = read_csv(path, sep = "\t") #I preprocessed with excel. Might need to fix this later.

        return mqoutput_df

    def _remove_decoys(self, mqoutput_df):

        mqoutput_df = mqoutput_df[~mqoutput_df["Protein IDs"].str.startswith("REV")]
        mqoutput_df = mqoutput_df[~mqoutput_df["Protein IDs"].str.startswith("CON")]

        return mqoutput_df

    def get_major_proteins(self, mqoutput_df):
        
        #setting this is just an object thing, all values should be unique...
        majority_proteins_set = set(mqoutput_df["Protein IDs"].to_list())

        return majority_proteins_set

    def get_protein_groups_dictionary(self, mqoutput_df):
        '''Quicker if a little risky, assume that 
        major proteins only appear once. (seems true on test 
        data). Iterate over rows to get protein group dictionary'''

        protein_groups = dict()
        for _, row in mqoutput_df.iterrows():
            protein_groups[row["Majority protein IDs"]] = set(row["Protein IDs"].split(";"))


        return protein_groups