class ProcessedPSMs():
    """
    A class that bundles a psm table in a pandas dataframe
    with some useful methods.

    Attributes:
    -----------
    df : pandas.DataFrame
        A pandas dataframe containing the psms list.

    Methods:
    --------
    get_proteins()
        Returns a list of proteins contained in the dataframe
        attribute.
    get_peptides()
        Returns a list of peptides contained in the dataframe
        attribute.
    get_shape()
        Returns the shape of the dataframe attribute in a 
        length 2 tuple.

    """

    def __init__(self, df):
        self.df = df

    def get_proteins(self):
        return set(self.df["protein id"])

    def get_peptides(self):
        return set(self.df["sequence"])

    def get_shape(self):
        return self.df.shape
