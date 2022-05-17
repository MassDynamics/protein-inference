from pandas import Series, DataFrame
from networkx import Graph
from multiprocessing import cpu_count, Pool
from operator import itemgetter

class ProteinMerger():
    """
    A class responsible for identifying indistinguishable proteins
    inside a problem network and mergine those nodes.

    Methods:
    --------
    run(pn):
        Merges proteins with the the same peptide connections and 
        returns an updated problem network.
    _list_to_string(l):
        Utility for converting a list to a string seperated by
        semi-colons. 
    _string_to_list(l):
        Utility for a string seperated by
        semi-colons to a list.
    get_mapping(pn):
        Identifies proteins with same peptide neighbors, returning
        a dataframe of peptide neighbour list, protein list pairs. 
    get_named_proteins(pn):
        Chooses proteins from the protein lists with the same
        peptide neighbours to keep (by sort order). Returns
        a dataframe indicating these.

    """

    def run(self, pn):

        mapping_df = self.get_named_proteins(pn)

        for _, row in mapping_df.iterrows():
            # Use protein 0
            named_protein = row["named"]
            pn.network.nodes[named_protein]["indistinguishable"] = sorted(row["protein"])
            # delete others
            pn.network = Graph(pn.network)
            to_remove = row["protein"]
            pn.network.remove_nodes_from(to_remove)
                

        return pn

    def _list_to_string(self, l):
        return ";".join(l)

    def _string_to_list(self, s):
        return s.split(";")

    def get_mapping(self, pn):

        s1 = Series(pn.get_proteins())
        s2 = s1.apply(lambda x: list(pn.network.neighbors(x)))
        s2 = s2.apply(lambda x: self._list_to_string(sorted(x)))
        df = DataFrame([s1, s2]).T
        df.columns = ["protein", "peptides"]
        df = df.groupby("peptides").agg({"protein": list}).reset_index()
        df.peptides = df.peptides.apply(lambda x: self._string_to_list(x))
        df.protein = df.protein.apply(lambda x: sorted(x))
        return df

    def get_named_proteins(self, pn):

        mapping_df = self.get_mapping(pn)
        named_proteins = []
        for _, row in mapping_df.iterrows():
            score_dict = dict(zip(row["protein"], [pn.network.nodes[x]["score"] for x in row["protein"]]))
            best_scoring_protein = max(score_dict.items(), key = itemgetter(1))[0]
            #deal with same score proteins! (id's should be unique even if scores aren't)
            best_scoring_proteins = [k for k,v in score_dict.items() if v == score_dict[best_scoring_protein]]
            named_proteins.append(sorted(best_scoring_proteins)[0])

            row["protein"].remove(best_scoring_protein)

        mapping_df["named"] = named_proteins
        return mapping_df
