from scipy.stats import chi2
from numpy import log
from multiprocessing import cpu_count, Pool

#unfinished class which requires peptide level p-values are present 
#these aren't produced by percolator output but if you were using 
#psms from elsewhere then you might want to use this. 
#def fishers_method(self, pn):
#
#    for protein in pn.get_proteins():
#        k = 0 
#        stat =0
#        for peptide in pn.network.neighbors(protein):
#            k = k + 1
#            #this should be p-value, are we even loading this data?
#            stat = stat - 2*log(pn.network.edges[(protein,peptide)]["q_value"],2)
#        pn.network.nodes[protein]["score"] = 1 - chi2.cdf(stat,k)
#    
#    return pn

class PEPProductScorer():
    """A class that includes one function, run, which takes a 
    problem network and scores peptides use the PEP Product
    method."""

    def run(self, pn):
        """
        Scores proteins using the PEP Product Method. 

        Parameters
        ----------
        pn : ProblemNetwork
            An instance of the problem network object

        Returns
        -------
        pn : ProblemNetwork
            An instance of the problem network object where
            protein nodes have a score attribute.
        """
        for protein in pn.get_proteins():
            score = 1
            for peptide in pn.network.neighbors(protein):
                score = score*pn.network.edges[(protein,peptide)]["PEP"]
            pn.network.nodes[protein]["score"] = score

        return pn

class BestPeptideScorer():
    """A class that includes one function, run, which takes a 
    problem network and scores peptides use the best peptide
    method."""

    def run(self, pn):
        """
        Scores proteins using the PEP Product Method. 

        Parameters
        ----------
        pn : ProblemNetwork
            An instance of the problem network object

        Returns
        -------
        pn : ProblemNetwork
            An instance of the problem network object where
            protein nodes have a score attribute.
        """
        for protein in pn.get_proteins():
            scores = []
            for peptide in pn.network.neighbors(protein):
                scores.append(pn.network.edges[(protein,peptide)]["score"])
            
            #every protein will have at least one score
            pn.network.nodes[protein]["score"] = sorted(scores)[-1]

        return pn

class BestTwoPeptideScorer:
    """A class that includes one function, run, which takes a 
    problem network and scores peptides use the best peptide
    method."""

    def run(self, pn):
        """
        Scores proteins using the Best Two Peptides method. 

        Parameters
        ----------
        pn : ProblemNetwork
            An instance of the problem network object

        Returns
        -------
        pn : ProblemNetwork
            An instance of the problem network object where
            protein nodes have a score attribute.
        """

        for protein in pn.get_proteins():
            scores = []
            for peptide in pn.network.neighbors(protein):
                scores.append(pn.network.edges[(protein,peptide)]["score"])
            
            if len(scores) > 1:
                pn.network.nodes[protein]["score"] = sorted(scores)[-1] + sorted(scores)[-2]
            
            else:
                pn.network.nodes[protein]["score"] = 0
        return pn