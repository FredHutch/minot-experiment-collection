
import boto3
import gzip
import io
import json
import logging
import numpy as np
import os
import pandas as pd
import shutil
import subprocess
import sys
import traceback

from scipy.stats.mstats import gmean


def repack_hdf5(fp, filter_string="GZIP=7"):
    """Repack an HDF5 file."""
    # Use a temporary file to store the repacked HDF5 file
    temp_fp = fp + '.repacked.hdf5'
    
    # Format the system call command
    commands = [
        "h5repack",
        "-i", fp,
        "-o", temp_fp,
        "-f", filter_string
    ]
    
    logging.info("Running command: {}".format(" ".join(commands)))

    # Run the command
    p = subprocess.Popen(
        commands,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )
    stdout, stderr = p.communicate()
    exitcode = p.wait()

    if stdout:
        logging.info("Standard output of subprocess:")
        for line in stdout.decode("utf-8").split('\n'):
            logging.info(line)
    if stderr:
        logging.info("Standard error of subprocess:")
        for line in stderr.split('\n'):
            logging.info(line)

    # Check the exit code
    if exitcode != 0:
        assert exitcode == 0, "Exit code {}".format(exitcode)

    # Copy the repacked HDF5 back to the original file path
    logging.info("Copying {} to {}".format(temp_fp, fp))
    shutil.copyfile(temp_fp, fp)


def format_eggnog_ko_df(df):
    """Make a table with just genes and KOs."""

    return pd.DataFrame([
        {
            "gene": r["#query_name"],
            "ko": ko
        }
        for _, r in df.loc[
            df["KEGG_KOs"].apply(
                lambda v: isinstance(v, str) and len(v) > 0
                )
            ].iterrows()
        for ko in r["KEGG_KOs"].split(",")
    ])


def format_eggnog_go_df(df):
    """Make a table with just genes and GOs."""
    return pd.DataFrame([
        {
            "gene": r["#query_name"],
            "go": go
        }
        for _, r in df.loc[
            df["GO_terms"].apply(
                lambda v: isinstance(v, str) and len(v) > 0
            )
        ].iterrows()
        for go in str(r["GO_terms"]).split(",")
    ]).dropna()


def format_eggnog_cluster_df(df):
    """Make a table with just genes and eggNOG clusters."""
    return df.loc[:, ["#query_name", "seed_eggNOG_ortholog"]].rename(columns={
        "seed_eggNOG_ortholog": "eggnog_cluster",
        "#query_name": "gene"
    }).dropna()


def add_table_to_store(
    metadata_table_fp,
    store,
    filter_function_dict,
    sep="\t",
    header="infer",
    names=None,
    data_columns=None,
    comment=None
):
    """
    Add a table to the store.

    In `filter_function_dict`, each key is the name of the table, 
    and the value is a function to transform the table.

    This is somewhat convoluted, but it is one way to maximize code reuse.
    
    """

    # Read in the table with Pandas
    logging.info("Reading in {}".format(metadata_table_fp))
    df = pd.read_table(metadata_table_fp, sep=sep, header=header, names=names, comment=comment)

    # Replace the NaN values with "none" to prevent errors writing to HDF5
    df.fillna("none", inplace=True)

    # Use each of the filter functions to write out a table to the store
    for table_name, filter_function in filter_function_dict.items():
        logging.info("Applying filter function for {}".format(table_name))
        filtered_df = filter_function(df)

        logging.info("Writing a table with {:,} rows and {:,} columns to HDF5".format(
            filtered_df.shape[0],
            filtered_df.shape[1]
        ))

        if data_columns is not None:
            filtered_data_columns = [ix for ix in data_columns if ix in filtered_df.columns.values]
            logging.info("Using data columns: {}".format(", ".join(filtered_data_columns)))
        else:
            filtered_data_columns = None

        # Do not add any NaNs
        assert filtered_df.shape[0] == filtered_df.dropna().shape[0], "NaNs in DataFrame"

        filtered_df.to_hdf(
            store,
            table_name,
            format="table",
            data_columns=filtered_data_columns
        )


def add_cags_to_store(cags_json, store):
    """Add a set of CAGs to the HDF5 file as a table."""
    cags = read_json(cags_json)
    
    m = "CAGs must be formatted as a dict of lists"
    assert isinstance(cags, dict), m
    assert all([isinstance(v, list) for v in cags.values()]), m
    
    cags_df = pd.DataFrame([
        {
            "cag": cag_id,
            "gene": gene_id
        }
        for cag_id, gene_id_list in cags.items()
        for gene_id in gene_id_list
    ])
    
    assert cags_df.shape[0] > 0, "No CAGs were detected"

    cags_df.to_hdf(store, "cags", format="table", data_columns=["cag", "gene"])

    return cags


def add_abundance_to_store(sample_name, sample_abundance_json_fp, store, cags,
                           results_key="results", abundance_key="depth", other_keys=["length", "coverage", "nreads"], gene_id_key="id"):
    """Add the abundance data from a single abundance JSON to the store."""
    
    # Get the JSON for this particular sample
    sample_dat = read_json(sample_abundance_json_fp)

    # Make sure that the key for the results is in this file
    assert results_key in sample_dat

    # Subset down to the list of results
    sample_dat = sample_dat[results_key]
    assert isinstance(sample_dat, list)

    # Make sure that every element in the list has the indicated keys
    for d in sample_dat:
        if other_keys is not None:
            for k in other_keys:
                assert k in d
        assert abundance_key in d
        assert gene_id_key in d

    # Format as a DataFrame
    sample_dat = pd.DataFrame(
        [
            {
                gene_id_key: d[gene_id_key],
                abundance_key: float(d[abundance_key]),
                **{
                    k: d[k]
                    for k in other_keys
                }
            }
            for d in sample_dat
        ]
    )

    # Add the sample name
    sample_dat["sample"] = sample_name

    logging.info("Sample {} contains {} genes".format(sample_name, sample_dat.shape[0]))

    # If the abundance isn't a CLR, calculate the CLR

    if abundance_key != "clr":
        # Make sure there are all positive values before trying to calculate the CLR
        if sample_dat[abundance_key] <= 0:
            logging.info("Cannot calculate the CLR with values <= 0")
        else:
            logging.info("Calculating the CLR")

            # Get the geometric mean of the sample
            sample_gmean = gmean(sample_dat[abundance_key])

            # Calculate the CLR
            sample_dat["clr"] = (sample_dat[abundance_key] / sample_gmean).apply(np.log10)

    # Write to the HDF5
    logging.info("Writing {} to HDF5".format(sample_name))
    sample_dat.to_hdf(
        store,
        "abundance",
        format="table",
        data_columns=[gene_id_key, "sample"],
        append=True
    )

    # If the CAGs are provided, make a summary of their abundance
    if cags is not None:
        logging.info("Calculating CAG abundances")
        abund_dict = sample_dat.set_index(gene_id_key)[abundance_key].to_dict()

        # Calculate the mean abundance of each CAG in this sample
        cag_df = pd.DataFrame([
            {
                "cag_id": cag_id,
                abundance_key: np.mean([
                    abund_dict.get(gene_id, 0)
                    for gene_id in gene_id_list
                ])
            }
            for cag_id, gene_id_list in cags.items()
        ])

        # Remove the CAGs that were not detected at all
        cag_df = cag_df.loc[cag_df[abundance_key] > 0]

        # Add the sample name
        cag_df["sample"] = sample_name

        # Calculate the CLR
        if abundance_key != "clr":
            cag_df["clr"] = cag_df[abundance_key].apply(
                lambda v: np.log10(v / sample_gmean)
            )

        logging.info("Writing out the abundance for {:,} CAGs".format(
            cag_df.shape[0]
        ))
        cag_df.to_hdf(
            store,
            "cag_abundance",
            format="table",
            data_columns=["cag_id", "sample"],
            append=True
        )


    logging.info("Done reading in abundance for {}".format(sample_name))


def exit_and_clean_up(temp_folder):
    """Log the error messages and delete the temporary folder."""
    # Capture the traceback
    logging.info("There was an unexpected failure")
    exc_type, exc_value, exc_traceback = sys.exc_info()
    for line in traceback.format_tb(exc_traceback):
        logging.info(line)

    # Delete any files that were created for this sample
    logging.info("Removing temporary folder: " + temp_folder)
    shutil.rmtree(temp_folder)

    # Exit
    logging.info("Exit type: {}".format(exc_type))
    logging.info("Exit code: {}".format(exc_value))
    sys.exit(exc_value)


def read_json(fp):
    assert fp.endswith((".json", ".json.gz"))
    logging.info("Reading in " + fp)
    if fp.startswith("s3://"):
        # Parse the S3 bucket and key
        bucket_name, key_name = fp[5:].split("/", 1)

        # Connect to the S3 boto3 client
        s3 = boto3.client('s3')

        # Download the object
        retr = s3.get_object(Bucket=bucket_name, Key=key_name)

        if fp.endswith(".gz"):
            # Parse GZIP
            bytestream = io.BytesIO(retr['Body'].read())
            got_text = gzip.GzipFile(
                None, 'rb', fileobj=bytestream).read().decode('utf-8')
        else:
            # Read text
            got_text = retr['Body'].read().decode('utf-8')

        # Parse the JSON
        dat = json.loads(got_text)

    else:
        assert os.path.exists(fp)

        if fp.endswith(".gz"):
            dat = json.load(gzip.open(fp, "rt"))
        else:
            dat = json.load(open(fp, "rt"))

    # Make sure that the sample sheet is a dictionary
    assert isinstance(dat, dict)

    return dat
