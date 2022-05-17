from pandas import concat
import matplotlib.pyplot as plt
#import plotly_express as px


class FalseDiscoveryRateCalculator():
    """
    Packages related function for estimating false discovery
    rates and q - values. 

    """
    def FDR(self, threshold, target, decoy):
        '''
        
        Calculates the false discovery rate at a given 
        threshold based on the score distributions in 
        the target and decoy pandas tables. '
        
        Parameters
        ----------
        threshold : float
            A number between 0 and 1.
        target : pandas.DataFrame
            A target protein table produced by ProteinInferenceRunner().get_output()
        decoy : pandas.DataFrame
            A decoy protein table produced by ProteinInferenceRunner().get_output()
        
        '''

        target_positives = sum(target.score >= threshold)
        decoy_positives = sum(decoy.score >= threshold)
        false_discovery_rate = decoy_positives / target_positives

        return false_discovery_rate

    def tag_FDR(self, target, decoy, entrapment=False):
        '''
        Creates an FDR column in the target table. If
        entrapment is set to true it changes the column
        name (used in benchmarking.)

        Parameters
        ----------
        target : pandas.DataFrame
            A target protein table produced by ProteinInferenceRunner().get_output()
        decoy : pandas.DataFrame
            A decoy protein table produced by ProteinInferenceRunner().get_output()
        entrapment : bool
            Indicates whether or not to prefix the FDR column with "entrapment"

        '''

        if not entrapment:
            target.loc[:, "FDR"] = target.score.apply(
                lambda x: self.FDR(x, target, decoy))
        else:
            target.loc[:, "entrapmentFDR"] = target.score.apply(
                lambda x: self.FDR(x, target, decoy))

        return target

    def tag_q_value(self, target, decoy, entrapment=False):
        '''
        Create a q-value column for the target table
        based on the score distributions in the
        target and decoy pandas tables.

        Parameters
        ----------
        target : pandas.DataFrame
            A target protein table produced by ProteinInferenceRunner().get_output()
        decoy : pandas.DataFrame
            A decoy protein table produced by ProteinInferenceRunner().get_output()
        entrapment : bool
            Indicates whether or not to prefix the FDR column with "entrapment"

        '''
        if "FDR" not in target:
            target = self.tag_FDR(target, decoy, entrapment)

        if "entrapmentFDR" not in target and entrapment:
            target = self.tag_FDR(target, decoy, entrapment)

        target = target.sort_values("score")
        if not entrapment:
            q_values = [min(target.FDR[:i]) for i in range(1, len(target.FDR))]
        else:
            q_values = [min(target.entrapmentFDR[:i])
                        for i in range(1, len(target.FDR))]
        q_values = q_values + [q_values[-1]]
        
        if entrapment:
            target["q-value-entrapment"] = q_values
        else:
            target["q-value"] = q_values

        return target.sort_values("q-value")
