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

    def gene_abundance(self, genes=None, samples=None, metric=None):
        """
        
        Return a DataFrame with the abundance for a set of samples.

        If `metric` is None, return the abundance key that was used in the input. 
        Another option would be "clr".

        """

        # Default to returning all samples
        if samples is None:
            samples = self.all_samples
        else:
            assert isinstance(samples, list)
            for sample_id in samples:
                assert sample_id in self.all_samples, "{} is not a valid sample".format(sample_id)
        
        # Set the metric to return
        if metric is None:
            metric = self.abund_id_key

        # Store the data in a dict, keyed by sample
        df = {}

        # Iterate over the samples
        for sample_id in samples:
            # Read the abundance
            df[sample_id] = self.sample_gene_abundance(
                sample_id,
                metric=metric
            )
            
            # Subset to the genes of interest
            if genes is not None:
                df[sample_id] = df[sample_id].reindex(genes)

        # Format as a DataFrame
        return pd.DataFrame(df)

    @lru_cache(maxsize=512)
    def sample_gene_abundance(self, sample_id, metric=None):
        """
        
        Return a DataFrame with the abundance for a single sample.

        If `metric` is None, return the abundance key that was used in the input. 
        Another option would be "clr".

        """

        # Set the metric to return
        if metric is None:
            metric = self.abund_id_key

        # Read the abundance
        abund = pd.read_hdf(self.exp_col_fp, "/abundance/" + sample_id)

        for k in [self.gene_id_key, metric]:
            assert k in abund.columns.values, "Column {} not found for {}".format(
                k, sample_id)

        return abund.set_index(self.gene_id_key)[metric]

    @lru_cache(maxsize=512)
    def sample_cag_abundance(self, sample_id, metric=None):
        """
        
        Return a DataFrame with the abundance of CAGs for a single sample.

        If `metric` is None, return the abundance key that was used in the input. 
        Another option would be "clr".

        """

        # Set the metric to return
        if metric is None:
            metric = self.abund_id_key

        # Read the abundance
        abund = pd.read_hdf(self.exp_col_fp, "/cag_abundance/" + sample_id)

        for k in ["cag_id", metric]:
            assert k in abund.columns.values, "Column {} not found for {}".format(
                k, sample_id)

        return abund.set_index("cag_id")[metric]

    @lru_cache(maxsize=1)
    def metadata(self):
        """Return the metadata table."""
        return pd.read_hdf(self.exp_col_fp, "metadata")

    @lru_cache(maxsize=2)
    def eggnog_annotation(self, annot_type="ko"):
        """Return the entire set of eggNOG annotations, 'ko' or 'cluster'."""

        assert annot_type in ["ko", "cluster"]

        if annot_type == "ko":
            table_name = "eggnog_ko"
            col_name = "ko"

        elif annot_type == "cluster":
            table_name = "eggnog_cluster"
            col_name = "eggnog_cluster"

        return pd.read_hdf(
            self.exp_col_fp, 
            table_name
        ).set_index("gene")[col_name]

    @lru_cache(maxsize=1)
    def taxonomic_annotation(self):
        """Return the entire set of taxonomic annotations."""

        return pd.read_hdf(
            self.exp_col_fp,
            "taxonomic_classification"
        ).set_index("gene")

    def cag_abundance(self, cags=None, samples=None, metric=None):
        """
        
        Return a DataFrame with the abundance of CAGs for a set of samples.

        If `metric` is None, return the abundance key that was used in the input. 
        Another option would be "clr".

        """

        if samples is None:
            samples = self.all_samples

        # Set the metric to return
        if metric is None:
            metric = self.abund_id_key

        # Store the data in a dict, keyed by sample
        df = {}

        # Iterate over the samples
        for sample_id in samples:
            # Read the abundance
            df[sample_id] = self.sample_cag_abundance(sample_id, metric=metric)

            if cags is not None:
                df[sample_id] = df[sample_id].reindex(cags)

        # Format as a DataFrame
        return pd.DataFrame(df)

    @lru_cache(maxsize=1)
    def cag_membership(self):
        """Return a dict with the genes in each CAG."""
        cags = pd.read_hdf(self.exp_col_fp, "cags")

        return {
            cag_id: cag_df["gene"].tolist()
            for cag_id, cag_df in cags.groupby("cag")
        }

    @lru_cache(maxsize=128)
    def contigs_with_gene(self, gene_id):
        """Get the list of contigs that contain a given gene."""
        return pd.read_hdf(
            self.exp_col_fp, 
            "gene_positions", 
            where="cluster == '{}'".format(gene_id)
        )["seqname"].tolist()

    @lru_cache(maxsize=128)
    def contig_df(self, contig_id):
        """Get the summary of the structure of a contig."""
        return pd.read_hdf(
            self.exp_col_fp,
            "gene_positions",
            where="seqname == '{}'".format(contig_id)
        )
