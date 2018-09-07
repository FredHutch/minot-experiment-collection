#!/bin/bash

./make-experiment-collection.py \
    --output-hdf5 test-experiment-collection.hdf5 \
    --output-logs test-experiment-collection.log \
    --abundance-sample-sheet tests/data/small_demonstration_experiment_2018.sample_sheet.local.json \
    --metadata-table tests/data/metadata.csv \
    --metadata-field-sep "," \
    --taxonomic-classification-tsv tests/data/small_demonstration_experiment_2018.nr.tax.20180717.diamond.tax.gz \
    --eggnog-mapper-tsv tests/data/small_demonstration_experiment_2018.eggnog-mapper-20180601.tsv.gz \
    --integrated-assembly tests/data/small_demonstration_experiment_2018.hdf5 \
    --cags-json tests/data/small_demonstration_experiment_2018_2_samples_clr_0.05.cags.json.gz \
    --temp-folder ./
