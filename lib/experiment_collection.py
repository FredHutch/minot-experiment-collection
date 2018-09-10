"""Class to help read data from the experiment collection."""

from collections import defaultdict
import os
import pandas as pd

from functools import lru_cache

class ExperimentCollection:

    def __init__(self, exp_col_fp, gene_id_key="id", abund_id_key="depth"):
        """Pass in the filepath for the experiment collection."""

        # Save the filepath
        self.exp_col_fp = exp_col_fp

        # Make sure the file exists
        assert os.path.exists(self.exp_col_fp)

        # Set the default gene ID key
        self.gene_id_key = gene_id_key
        # Set the default abundance key
        self.abund_id_key = abund_id_key

        # Get the list of all samples that have abundance information
        with pd.HDFStore(self.exp_col_fp, mode="r") as store:
            self.all_samples = [
                sample_id.replace("/abundance/", "")
                for sample_id in store
                if sample_id.startswith("/abundance/")
            ]

    @lru_cache(maxsize=2)
    def gene_abundance(self, metric=None):
        """
        
        Return a DataFrame with the abundance for a set of samples.

        If `metric` is None, return the abundance key that was used in the input. 
        Another option would be "clr".

        """
        
        # Open a connection to the HDFStore
        store = pd.HDFStore(self.exp_col_fp, mode="r")

        # Set the metric to return
        if metric is None:
            metric = self.abund_id_key

        # Store the data in a dict, keyed by sample
        df = {}

        # Iterate over the samples
        for sample_id in self.all_samples:
            # Read the abundance
            abund = pd.read_hdf(store, "/abundance/" + sample_id)

            for k in [self.gene_id_key, metric]:
                assert k in abund.columns.values, "Column {} not found for {}".format(k, sample_id)

            # Add to the dict
            df[sample_id] = abund.set_index(self.gene_id_key)[metric]

        # Close the connection to the HDF5 file
        store.close()

        # Format as a DataFrame
        return pd.DataFrame(df)

    @lru_cache(maxsize=1)
    def metadata(self):
        """Return the metadata table."""
        return pd.read_hdf(self.exp_col_fp, "metadata")

    @lru_cache(maxsize=2)
    def eggnog_annotation(self, annot_type="ko"):
        """Return the entire set of eggNOG annotations, 'ko' or 'cluster'."""

        assert annot_type in ["ko", "cluster"]

        if annot_type == "ko":

            return pd.read_hdf(self.exp_col_fp, "eggnog_ko")

        elif annot_type == "cluster":

            return pd.read_hdf(self.exp_col_fp, "eggnog_cluster")

    @lru_cache(maxsize=1)
    def taxonomic_annotation(self):
        """Return the entire set of taxonomic annotations."""

        return pd.read_hdf(self.exp_col_fp, "taxonomic_classification")

    @lru_cache(maxsize=2)
    def cag_abundance(self, metric=None):
        """
        
        Return a DataFrame with the abundance of CAGs for a set of samples.

        If `metric` is None, return the abundance key that was used in the input. 
        Another option would be "clr".

        """

        # Open a connection to the HDFStore
        store = pd.HDFStore(self.exp_col_fp, mode="r")

        # Set the metric to return
        if metric is None:
            metric = self.abund_id_key

        # Store the data in a dict, keyed by sample
        df = {}

        # Iterate over the samples
        for sample_id in self.all_samples:
            # Read the abundance
            abund = pd.read_hdf(store, "/cag_abundance/" + sample_id)
            
            for k in ["cag_id", metric]:
                assert k in abund.columns.values, "Column {} not found for {}".format(
                    k, sample_id)

            # Add to the dict
            df[sample_id] = abund.set_index("cag_id")[metric]

        # Close the connection to the HDF5 file
        store.close()

        # Format as a DataFrame
        return pd.DataFrame(df)
