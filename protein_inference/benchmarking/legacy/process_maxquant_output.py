import os
from os.path import join
import pandas as pd 
import pickle
from pandas import read_csv
from protein_inference.processing.psms_preprocessor import PSMsPreprocessor
from protein_inference.processing.psms_network_splitter import PSMsNetworkSplitter
from protein_inference.processing.psms_network_generator import PSMsNetworkGenerator
from protein_inference.table_maker import TableMaker

class ProcessMaxQuant():

    def get_peptide_unique_mapping(self, mq_peptide_table):
        df =  mq_peptide_table[["Sequence", 'Unique (Proteins)']]
        return self.slow_multi_mapping(df,"Sequence",'Unique (Proteins)')

    def get_peptide_razor_mapping(self, mq_peptide_table):
        df =  mq_peptide_table[["Sequence", 'Unique (Groups)']]
        return self.slow_multi_mapping(df,"Sequence",'Unique (Groups)')

    def get_peptide_group_mapping(self,mq_peptide_table):
        mask = mq_peptide_table["Protein group IDs"].str.split(";").apply(lambda x: len(x))!=1
        mq_peptide_table.loc[mask, "Protein group IDs"] = "multiple"
        return self.slow_multi_mapping(mq_peptide_table,"Sequence","Protein group IDs")

    def get_protein_score_mapping(self,mq_protein_table):
        df = mq_protein_table
        df = df.assign(protein_id=df['Majority protein IDs'].str.split(';')).explode('protein_id')
        return self.slow_multi_mapping(df,"protein_id","Score")

    def get_protein_group_mapping(self, mq_protein_table):
        df = mq_protein_table[["Protein IDs"]].reset_index()
        df = df.assign(protein_id=df['Protein IDs'].str.split(';')).explode('protein_id')
        return self.slow_multi_mapping(df,"protein_id","index")

    def slow_multi_mapping(self, df, key, val):
        '''A utility for creating 1 to many mappings from pandas dataframes.
        Probably very slow.'''
        def get_mapping(i):
            mask = df[key] == i
            value = df[val][mask]
            if len(value.values) > 0:
                return value.values[0]
            else:
                return "not found"
        return get_mapping

    def make_psms_standard(self, mq_psms_table, decoy = False):

        #talk to Peppe about this
        #mq_psms_table.Sequence = mq_psms_table["Modified sequence"]

        #remove decoy proteins
        mq_psms_table = mq_psms_table[mq_psms_table.Proteins.str.contains("REV_") == decoy]
        #remove reversed peptides
        mq_psms_table = mq_psms_table[mq_psms_table.Reverse != "+"]
        #remove other columns
        mq_psms_table = mq_psms_table[["Sequence","Proteins","PEP","Score"]]
        #set q_value to missing
        mq_psms_table["q_value"] = -1
        #change multiple protein matches seperator
        mq_psms_table.Proteins = mq_psms_table.Proteins.str.replace(';',',')

        return mq_psms_table

    def run(self, mq_psms_table, mq_peptide_table, mq_protein_table):

        #create network
        mq_psms_table = self.make_psms_standard(mq_psms_table)
        psms = PSMsPreprocessor(mq_psms_table).get_processed_psms()
        network = PSMsNetworkGenerator(psms).generate_network()

        #get score and group mappings
        pep_unique_map = self.get_peptide_unique_mapping(mq_peptide_table)
        pep_razor_map =  self.get_peptide_razor_mapping(mq_peptide_table)
        pep_group_map =  self.get_peptide_group_mapping(mq_peptide_table)
        pro_group_map =  self.get_protein_group_mapping(mq_protein_table)
        pro_score_map =  self.get_protein_score_mapping(mq_protein_table)

        #need this to get group name from group index
        majority_from_group = mq_protein_table["Majority protein IDs"].str.split(';').apply(lambda x: sorted(x)[0])
        #indistinguishable_from_group = mq_protein_table["Majority protein IDs"].str.split(';').apply(lambda x: sorted(x)[1:])

        for protein in psms.get_proteins():
            network.nodes[protein]["major"] = majority_from_group[pro_group_map(protein)]
            network.nodes[protein]["unique_evidence"] = False #until proven otherwise
            #only record indistinguishable for major
            #if protein == majority_from_group[pro_group_map(protein)]:
            #    network.nodes[protein]["indistinguishable"] = indistinguishable_from_group[int(pro_group_map(protein))]
            
            score = pro_score_map(protein)
            if type(score) is not str:
                network.nodes[protein]["score"] = score
            else:
                network.nodes[protein]["score"] = 0

        for peptide in mq_peptide_table.Sequence.unique():
            if peptide in network.nodes():
                #set uniqueness
                network.nodes[peptide]["unique"] = pep_unique_map(peptide) == "yes"
                #set razorness
                network.nodes[peptide]["razor"] = pep_razor_map(peptide) == "yes"
                
                #set group
                if pep_group_map(peptide) == 'multiple':
                    network.nodes[peptide]["allocated"] = "Multiple"
                else:
                    network.nodes[peptide]["allocated"] = majority_from_group[int(pep_group_map(peptide))]

                #if unique, set unique evidenced to true in representation:
                if pep_unique_map(peptide) == "yes":
                    protein = list(network.neighbors(peptide))[0]
                    network.nodes[protein]["unique_evidence"] = True

        # set protein unique evidenced:
        problem_networks = PSMsNetworkSplitter(network).split_networks()

        # create protein table
        protein_table = TableMaker().get_system_protein_table(problem_networks)
        protein_table = protein_table.sort_values(["score","protein_id"], ascending = [False,True]).reset_index(drop = True)

        # label protein q-value (call FDR for now)
        qvalue_dict = self.get_q_value_dict(mq_protein_table)
        protein_table["q-value"] =  protein_table.protein_id.apply(lambda x: self.get_q_value(qvalue_dict, x))

        # create peptide table
        peptide_table = TableMaker().get_system_peptide_table(problem_networks)

        return problem_networks, protein_table, peptide_table

    def get_q_value_dict(self, mq_protein_table):
        df = mq_protein_table[['Majority protein IDs',"Q-value"]]
        df = df.assign(id = df['Majority protein IDs'].str.split(';')).explode("id")
        return dict(zip(df.id, df["Q-value"]))

    def get_q_value(self, qvalue_dict, id):
        if id in qvalue_dict.keys():
            return qvalue_dict[id]
        return 1 #worst possible q-value for proteins not assigned to "majority proteins"


class ReformatMaxQuant():

    def get_mq_inference(self, mq_home):

        # read data
        mq_psms_table = pd.read_csv(join(mq_home, "msms.txt"), sep = "\t").set_index("Peptide ID")
        mq_peptide_input = pd.read_csv(join(mq_home, "peptides.txt"), sep = "\t")
        mq_protein_input = pd.read_csv(join(mq_home, "ProteinGroups.txt"), sep = "\t")

        return mq_psms_table, mq_peptide_input, mq_protein_input

    def run(self, mq_home):

        # read data
        mq_psms_table, mq_peptide_input, mq_protein_input = self.get_mq_inference(mq_home)

        # process data
        mq_networks, mq_protein_table, mq_peptide_table = ProcessMaxQuant().run(mq_psms_table, mq_peptide_input, mq_protein_input)

        # write data
        self.write_mq_reformatted(mq_home, mq_networks, mq_protein_table, mq_peptide_table)

        return 1


    def write_mq_reformatted(self, mq_home, mq_networks, mq_protein_table, mq_peptide_table):
        # write data
        if not os.path.exists(os.path.join(mq_home, "reprisal_format")):
            os.makedirs(os.path.join(mq_home, "reprisal_format"))
        mq_protein_table.to_csv(os.path.join(mq_home, "reprisal_format","protein_table.csv"),index=False)
        mq_peptide_table.to_csv(os.path.join(mq_home, "reprisal_format","peptide_table.csv"),index=False)
        pickle.dump(mq_networks, open(os.path.join(mq_home, "reprisal_format","mq_networks.p"), "wb"))
        return 1