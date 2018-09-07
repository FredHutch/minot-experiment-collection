# Metagenomic (WGS) Experiment Collection

Author: Samuel S. Minot, Ph.D.

Make a single HDF5 with the complete set of data from a metagenomic (WGS) microbiome experiment

usage: make-experiment-collection.py [-h] --output-hdf5 OUTPUT_HDF5
                                     --output-logs OUTPUT_LOGS
                                     [--abundance-sample-sheet ABUNDANCE_SAMPLE_SHEET]
                                     [--cags-json CAGS_JSON]
                                     [--metadata-table METADATA_TABLE]
                                     [--metadata-field-sep METADATA_FIELD_SEP]
                                     [--taxonomic-classification-tsv TAXONOMIC_CLASSIFICATION_TSV]
                                     [--eggnog-mapper-tsv EGGNOG_MAPPER_TSV]
                                     [--integrated-assembly INTEGRATED_ASSEMBLY]
                                     [--temp-folder TEMP_FOLDER]

Collect all available information about a microbiome metagenomic WGS
experiment.

optional arguments:
  -h, --help            show this help message and exit
  --output-hdf5 OUTPUT_HDF5
                        Location for output HDF5 file.
  --output-logs OUTPUT_LOGS
                        Location for output logs from running this script.
  --abundance-sample-sheet ABUNDANCE_SAMPLE_SHEET
                        Location for sample sheet listing abundance files
                        (e.g. FAMLI output) [.json[.gz]].
  --cags-json CAGS_JSON
                        Location of JSON describing CAG membership for each
                        gene.
  --metadata-table METADATA_TABLE
                        Location of metadata table.
  --metadata-field-sep METADATA_FIELD_SEP
                        Field separator for metadata table.
  --taxonomic-classification-tsv TAXONOMIC_CLASSIFICATION_TSV
                        Location of TSV with output of DIAMOND taxonomic
                        assignment for each gene.
  --eggnog-mapper-tsv EGGNOG_MAPPER_TSV
                        Location of TSV with output of eggNOG mapper for each
                        gene.
  --integrated-assembly INTEGRATED_ASSEMBLY
                        Location of HDF5 with information on the integrated
                        assembly.
  --temp-folder TEMP_FOLDER
                        Folder for temporary files.