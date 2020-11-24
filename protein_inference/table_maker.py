from pandas import concat, DataFrame
from multiprocessing import cpu_count, Pool
from copy import deepcopy
import pandas as pd

class TableMaker():  
    """
    A class that groups functions required to generate output tables. 
    
    Methods:
    --------
    get_protein_table(pn)
        Processes a ProblemNetwork to retrieve the corresponding protein
        table
    get_protein_tables(pns)
        Processes a list of ProblemNetworks to retrieve the corresponding protein
        tables in a list.
    get_system_protein_table(pns):
        Processes a list of ProblemNetworks to retrieve the corresponding protein
        table.
    get_peptide_table(pn)
        Processes a ProblemNetwork to retrieve the corresponding peptide
        table
    get_peptide_tables(pns)
        Processes a list of ProblemNetwork to retrieve the corresponding peptide
        tables in a list.
    get_system_peptide_table(pns):
        Processes a list of ProblemNetworks to retrieve the corresponding peptide
        table.
    _get_edge_list(pn):
        Retrieves a pandas dataframe encoding the edges from peptide to protein 
        for each edge in a ProblemNetwork.
    find_molecule(pns, molecule):
        Find the molecule (a string corresponding to a protein id or modified
        peptide sequence) in the list of ProblemNetworks and returns the 
        corresponding network. 
        
    """

    def get_protein_table(self, pn):
        # source the data as a list
        df = self._get_edge_list(pn)

        # perform protein level aggregation
        agg_dict = self._get_agg_dict(pn)
        protein_df = df.groupby("protein_id").agg(agg_dict)

        # calculate new columns
        protein_df["total_peptides"] = df.groupby("protein_id").size()

        #if we know unique
        if pn.get_node_attribute_dict("unique"):
            protein_df["non_unique"] = protein_df.total_peptides - protein_df.unique 
            protein_df.non_unique = protein_df.non_unique.astype("int32")
        
        for col in ["razor","non_unique","unique"]:
            if pn.get_node_attribute_dict(col):
                protein_df[col] = protein_df[col].astype("int32")

        if pn.get_node_attribute_dict("non_unique"):
            protein_df.non_unique = protein_df.non_unique.astype("int32")

        if pn.get_node_attribute_dict("non_unique"):
            protein_df.non_unique = protein_df.non_unique.astype("int32")

        #sort sequence modified
        protein_df.sequence_modified = protein_df.sequence_modified.apply(lambda x: sorted(x))

        protein_df = protein_df.reset_index()  # otherwise protein id isn't a column

        # if solved, add new scores:
        dict_score = pn.get_node_attribute_dict("score")
        if dict_score:
            protein_df["score"] = protein_df.protein_id.apply(
                lambda x: dict_score[x])

        # if solved, add subset proteins:
        dict_subset = pn.get_node_attribute_dict("major")
        if dict_subset:
            protein_df["Group"] = protein_df.protein_id.apply(
                lambda x: dict_subset[x])

        if pn.get_node_attribute_dict("indistinguishable"):
            protein_df = self.add_indistinguishable_col(pn, protein_df)
            protein_df.indistinguishable = protein_df.indistinguishable.apply(lambda x: sorted(x))
        
        if pn.get_node_attribute_dict("major"):
            protein_df = self.add_subset_col(pn, protein_df)
            protein_df.subset = protein_df.subset.apply(lambda x: sorted(x))        

        cols = ["protein_id", 
                "unique","non_unique", "razor", "total_peptides",
                "score", "q-value",
                "Group", "indistinguishable", "subset",
                "sequence_modified"]

        new_cols = []
        for col in cols:
            if col in protein_df.columns:
                new_cols.append(col)
        
        protein_df = protein_df.loc[:,new_cols]

        if "score" in protein_df.columns:
            return protein_df.sort_values("score", ascending=False)
        else:
            return protein_df

    def get_protein_tables(self, pns):

        p = Pool(cpu_count())
        protein_tables = p.map(self.get_protein_table, pns)

        return protein_tables

    def get_system_protein_table(self, pns):
        protein_table = concat(self.get_protein_tables(pns))
        protein_table = self.emulate_percolator_formatting(protein_table)
        return protein_table

    def get_peptide_table(self, pn):

        def label_score(protein, score_dict):
            return score_dict[protein]

        protein_table = self.get_protein_table(pn)
        score_dict = dict(zip(protein_table.protein_id, protein_table.score))

        df = self._get_edge_list(pn)
        df["protein_score"] = df.apply(lambda row: label_score(row["protein_id"],
                                                    score_dict), axis=1)

        df = df.sort_values(["sequence_modified", "protein_score"],
                            ascending=[True, False])

        # remove duplicate proteins
        df = df.groupby("sequence_modified").first().reset_index()
        
        df = df.rename(columns = {"q_value":"q-value","sequence_modified":"sequence"})
        cols = ['sequence', 'unique','razor', 'PEP', 'q-value', 'score', 'major','protein_score', 'unique_evidence']

        new_cols = []
        for col in cols:
            if col in df.columns:
                new_cols.append(col)
        
        df = df.loc[:,new_cols]
        return df

    def get_peptide_tables(self, pns):

        p = Pool(cpu_count())
        peptide_tables = p.map(self.get_peptide_table, pns)

        return peptide_tables

    def get_system_peptide_table(self, pns):
        peptide_table = concat(self.get_peptide_tables(pns))
        return peptide_table

    def _get_edge_list(self, pn):

        rows = []
        for u, v, d in pn.network.edges(data=True):
            node_1_data = pn.network.nodes[u]
            node_2_data = pn.network.nodes[v]
            row = dict()
            if node_1_data["protein"]:
                row["protein_id"] = u
                row["sequence_modified"] = v
            else:
                row["protein_id"] = v
                row["sequence_modified"] = u
            row.update(node_1_data)
            row.update(node_2_data)
            row.update(d)
            rows.append(row)

        df = DataFrame(rows)

        return df.drop("protein", axis=1)

    # does this method belong here? It's more of a general utility
    def find_molecule(self, pns, molecule):
        for pn in pns:
            if molecule in pn.get_proteins():
                return pn
            elif molecule in pn.get_peptides():
                return pn
            else:
                indist_dict = pn.get_node_attribute_dict("indistinguishable")
                if indist_dict:
                    for l in indist_dict.values():
                        if molecule in l:
                            return pn

    def _flip_dict(self, old_dict):
        '''https://www.geeksforgeeks.org/python-program-to-swap-keys-and-values-in-dictionary/'''
        new_dict = {} 
        for key, value in old_dict.items(): 
            if value in new_dict: 
                new_dict[value].append(key) 
            else: 
                new_dict[value]=[key] 
        return new_dict

    def add_indistinguishable_col(self, pn, table):

        indistinguishable_dict = pn.get_node_attribute_dict("indistinguishable")
        new_col = []
        for _, row in table.iterrows():
            new_col.append(indistinguishable_dict[row["protein_id"]])
        
        table["indistinguishable"] = new_col

        return table

    def add_subset_col(self, pn, table):
        
        subset_dict = pn.get_node_attribute_dict("major")
        subset_dict = self._flip_dict(subset_dict)
        new_col = []
        for _, row in table.iterrows():
            if row["protein_id"] == row["Group"]: #and row["protein_id"] in subset_dict.keys():
                new_col.append(subset_dict[row["protein_id"]])
            else:
                new_col.append([])
        
        table["subset"] = new_col

        return table

    def _get_agg_dict(self,pn):
        agg_dict = {}
        if pn.get_node_attribute_dict("razor"):
            agg_dict.update({"razor":sum})
        if pn.get_node_attribute_dict("unique"):
            agg_dict.update({"unique":sum})
        if pn.get_node_attribute_dict("score"):
            agg_dict.update({"score":sum})

        #always add sequence_modified
        agg_dict.update({"sequence_modified":list})

        return agg_dict

    def emulate_percolator_formatting(self, protein_table):

        #rename columns
        col_dict = {"protein_id":"ProteinId","sequence_modified":"peptideIds"}
        protein_table = protein_table.rename(columns = col_dict)
        #sort before assigning groupID
        protein_table = protein_table.sort_values("ProteinId")
        #format peptides as percolator does
        protein_table["peptideIds"] = protein_table.peptideIds.apply(lambda x: " ".join(x))
        #add group ids
        labels, _ = pd.factorize(protein_table.Group)
        protein_table["ProteinGroupId"] = labels

        return protein_table