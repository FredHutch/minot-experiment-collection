#!/usr/bin/env python3
"""Collect all available information about a microbiome metagenomic WGS experiment."""

import argparse
import boto3
import logging
import os
import pandas as pd
import shutil
import sys
import uuid
from lib.helpers import exit_and_clean_up
from lib.helpers import read_json
from lib.helpers import add_abundance_to_store
from lib.helpers import add_cags_to_store
from lib.helpers import add_table_to_store
from lib.helpers import format_eggnog_cluster_df
from lib.helpers import format_eggnog_ko_df
# from lib.helpers import format_eggnog_go_df
from lib.helpers import repack_hdf5


def make_experiment_collection(
    output_hdf5=None,
    output_logs=None,
    abundance_sample_sheet=None,
    cags_json=None,
    metadata_table=None,
    metadata_field_sep=None,
    taxonomic_classification_tsv=None,
    eggnog_mapper_tsv=None,
    integrated_assembly=None,
    temp_folder=None
):

    # Make sure the temporary folder exists
    assert os.path.exists(temp_folder)

    # Make sure that at least one of the pieces of data has been specified
    assert any([x is not None for x in [
        abundance_sample_sheet, cags_json, metadata_table, taxonomic_classification_tsv,
        eggnog_mapper_tsv, integrated_assembly
    ]]), "No input data has been specified"

    # Make a new temp folder
    temp_folder = os.path.join(temp_folder, str(uuid.uuid4())[:8])
    os.mkdir(temp_folder)

    # Set up logging
    log_fp = os.path.join(temp_folder, "log.txt")
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [collect-experiment] %(message)s'
    )
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Write to file
    fileHandler = logging.FileHandler(log_fp)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    # Also write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    # Open a connection to AWS S3
    s3 = boto3.resource('s3')

    # If an integrated assembly HDF5 file was specified, copy it down and add to it
    local_hdf5_fp = os.path.join(temp_folder, "experiment.hdf5")

    if integrated_assembly is not None:
        logging.info("Copying integrated assembly from {}".format(integrated_assembly))

        if integrated_assembly.startswith("s3://"):
            bucket, key = integrated_assembly[5:].split("/", 1)
            try:
                s3.meta.client.download_file(bucket, key, local_hdf5_fp)
            except:
                exit_and_clean_up(temp_folder)
        else:
            try:
                assert os.path.exists(integrated_assembly)
            except:
                exit_and_clean_up(temp_folder)

            try:
                shutil.copyfile(integrated_assembly, local_hdf5_fp)
            except:
                exit_and_clean_up(temp_folder)

    # Add to that previous HDF5 file, if it exists, otherwise start a new one
    store = pd.HDFStore(local_hdf5_fp, mode="a")

    if cags_json is not None:
        logging.info("Reading in the CAGs and adding to the collection")

        try:
            cags = add_cags_to_store(cags_json, store)
        except:
            exit_and_clean_up(temp_folder)
    else:
        cags = None

    # Read in the sample_sheet
    if abundance_sample_sheet is not None:
        logging.info("Reading in the sample sheet from " + abundance_sample_sheet)
        try:
            abundance_sample_sheet = read_json(abundance_sample_sheet)
        except:
            exit_and_clean_up(temp_folder)

        try:
            assert isinstance(abundance_sample_sheet, dict), "Sample sheet must be a dict"
        except:
            exit_and_clean_up(temp_folder)

        logging.info("Adding sample abundance data to the collection")

        for sample_name, sample_abundance_json_fp in abundance_sample_sheet.items():
            for k in [".", "-"]:
                sample_name = sample_name.replace(k, "_")

            logging.info("Adding {} from {}".format(sample_name, sample_abundance_json_fp))

            try:
                add_abundance_to_store(sample_name, sample_abundance_json_fp, store, cags)
            except:
                exit_and_clean_up(temp_folder)

    if metadata_table is not None:
        logging.info("Reading in the metadata table and adding to the collection")

        try:
            add_table_to_store(
                metadata_table,
                store,
                {"metadata": lambda df: df},
                sep=metadata_field_sep
            )
        except:
            exit_and_clean_up(temp_folder)

    if taxonomic_classification_tsv is not None:
        logging.info(
            "Reading in the taxonomic classification table and adding to the collection")

        try:
            add_table_to_store(
                taxonomic_classification_tsv,
                store,
                {
                    "taxonomic_classification": lambda df: df.loc[df["taxid"] != 0, ["gene", "taxid"]]
                },
                sep="\t",
                header=None,
                names=["gene", "taxid", "evalue"],
                data_columns=["gene"]
            )
        except:
            exit_and_clean_up(temp_folder)

    if eggnog_mapper_tsv is not None:
        logging.info(
            "Reading in the eggNOG mapper results and adding to the collection")

        try:
            add_table_to_store(
                eggnog_mapper_tsv,
                store,
                {
                    "eggnog_ko": format_eggnog_ko_df,
                    # "eggnog_go": format_eggnog_go_df,
                    "eggnog_cluster": format_eggnog_cluster_df,
                },
                sep="\t",
                header=3,
                data_columns=["gene", "ko", "go", "eggnog_cluster"]
            )
        except:
            exit_and_clean_up(temp_folder)

    # Close the database
    store.close()

    # Repack and compress the database
    try:
        repack_hdf5(local_hdf5_fp)
    except:
        exit_and_clean_up(temp_folder)

    # Copy the file to the output
    for local_fp, remote_fp in [(local_hdf5_fp, output_hdf5), (log_fp, output_logs)]:
        logging.info("Copying {} to {}".format(
            local_fp, remote_fp
        ))
        if remote_fp.startswith("s3://"):
            bucket, key = remote_fp[5:].split("/", 1)
            s3.upload_file(local_fp, bucket, key)
        else:
            shutil.copyfile(local_fp, remote_fp)

    logging.info("Removing temporary folder")
    shutil.rmtree(temp_folder)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Collect all available information about a microbiome metagenomic WGS experiment."""
    )

    parser.add_argument("--output-hdf5",
                        type=str,
                        required=True,
                        help="""Location for output HDF5 file.""")
    parser.add_argument("--output-logs",
                        type=str,
                        required=True,
                        help="""Location for output logs from running this script.""")
    parser.add_argument("--abundance-sample-sheet",
                        type=str,
                        help="""Location for sample sheet listing abundance files (e.g. FAMLI output) [.json[.gz]].""")
    parser.add_argument("--cags-json",
                        type=str,
                        help="""Location of JSON describing CAG membership for each gene.""")
    parser.add_argument("--metadata-table",
                        type=str,
                        help="""Location of metadata table.""")
    parser.add_argument("--metadata-field-sep",
                        type=str,
                        default=",",
                        help="""Field separator for metadata table.""")
    parser.add_argument("--taxonomic-classification-tsv",
                        type=str,
                        help="""Location of TSV with output of DIAMOND taxonomic assignment for each gene.""")
    parser.add_argument("--eggnog-mapper-tsv",
                        type=str,
                        help="""Location of TSV with output of eggNOG mapper for each gene.""")
    parser.add_argument("--integrated-assembly",
                        type=str,
                        help="""Location of HDF5 with information on the integrated assembly.""")
    parser.add_argument("--temp-folder",
                        type=str,
                        default="/scratch",
                        help="Folder for temporary files.")

    args = parser.parse_args(sys.argv[1:])

    make_experiment_collection(
        **args.__dict__
    )
