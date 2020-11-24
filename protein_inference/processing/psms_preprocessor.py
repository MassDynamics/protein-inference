from protein_inference.processing.processed_psms import ProcessedPSMs


class PSMsPreprocessor:
    """
    A class that collects methods for preprocessing a psms list
    contained in a pandas dataframe.

    Attributes:
    -----------
    PSMs : pandas.DataFrame
        a pandas dataframe containing the psms list.  
    partial : int
        The number of psms to retrieve if a partial
        processing is to be done. 0 implies complete
        processing. 
    decoy : bool
        0 if input needs to be processed to remove
        decoy proteins

    Methods:
    --------
    get_processed_psms(df, partial=0, decoy=0)
        Applies a series of preprocessing methods to the
        PSMs attribute and returns the preprocessed 
        data as a ProcessedPSMs object.
    
    _preprocess_q_values(df)
        Returns only the psms with q_value < 0.01
    
    _preprocess_columns(df)
        Returns the highest scoring/lowest pep/lowest q value
        associated with each psm
    
    _preprocess_duplicates(df)
        Splits rows that have duplicate proteins id's
        into rows in the dataframe.
    
    _preprocess_decoys(df)
        Removes rows with a protein id prefixed by decoy

    _preprocess_column_names(df)
        Rename columns where percolator naming is used
        to more generic column names

    _preprocess_no_q_values(df)
        If no q values are present in the psms list, 
        set all q values to negative one.

    Warning: _preprocess_no_q_values(df) allows use of potentially
    bad psms because information is available for filtering
    matches. 
    """

    def __init__(self, df, partial=0, decoy=0):
        self.df = df
        self.partial = partial
        self.decoy = decoy

    def get_processed_psms(self):

        df = self.df.copy(deep=True)
        if self.partial != 0:
            print("Treated data as partial")
            df = df.sort_values("percolator score",
                                ascending=True).head(self.partial)  
        df = self._preprocess_column_names(df)                
        df = self._preprocess_q_value(df)
        df = self._preprocess_columns(df)
        df = self._preprocess_duplicates(df)
        if not self.decoy:
            df = self._preprocess_decoys(df)

        if df.shape[0] == 0:
            print("Warning: This might be a decoy dataset that you forgot to tell the preprocesser.")

        return ProcessedPSMs(df)

    def _preprocess_q_value(self, df):
        return df[df["q_value"] <= 0.01]

    def _preprocess_columns(self, df):
        agg_dict = {"score": "max",
                    "PEP": "min",
                    "q_value":"min"}
        df = df.groupby(["sequence", "protein id"]).agg(agg_dict).reset_index()

        return df

    def _preprocess_duplicates(self, df):
        # explode the protein id column
        df = df.assign(protein_id=df['protein id'].str.split(',')
                    ).explode('protein_id')

        df["protein id"] = df.protein_id
        df.drop("protein_id", axis=1)

        return df

    def _preprocess_decoys(self, df):
        return df[~df["protein id"].str.startswith('decoy')]

    def _preprocess_column_names(self, df):
        '''Autorename columns from percolator
        or other conventions formats'''

        if "percolator PEP" in df.columns:
            df = df.rename(columns = {"percolator PEP":"PEP"})

        if "percolator score" in df.columns:
            df = df.rename(columns = {"percolator score":"score"})

        if "percolator q-value" in df.columns:
            df = df.rename(columns = {"percolator q-value":"q_value"})
        
        if "Sequence" in df.columns:
            df = df.rename(columns = {"Sequence":"sequence"})

        if "Proteins" in df.columns:
            df = df.rename(columns = {"Proteins":"protein id"})

        if "Score" in df.columns:
            df = df.rename(columns = {"Score":"score"})

        return df

    def _preprocess_no_q_values(self,df):

        if "q_value" not in df.columns:
            df["q_value"] = -1

        return df
