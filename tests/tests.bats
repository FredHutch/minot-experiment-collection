#!/usr/bin/env bats

@test "make-experiment-collection.py in the PATH" {
  v="$(make-experiment-collection.py -h 2>&1 || true )"
  [[ "$v" =~ "Collect all available information" ]]
}

@test "Make experiment collection" {
  make-experiment-collection.py \
    --output-hdf5 test-experiment-collection.hdf5 \
    --output-logs test-experiment-collection.log \
    --abundance-sample-sheet /usr/local/tests/data/small_demonstration_experiment_2018.sample_sheet.Docker.json \
    --metadata-table /usr/local/tests/data/metadata.csv \
    --metadata-field-sep "," \
    --taxonomic-classification-tsv /usr/local/tests/data/small_demonstration_experiment_2018.nr.tax.20180717.diamond.tax.gz \
    --eggnog-mapper-tsv /usr/local/tests/data/small_demonstration_experiment_2018.eggnog-mapper-20180601.tsv.gz \
    --integrated-assembly /usr/local/tests/data/small_demonstration_experiment_2018.hdf5 \
    --cags-json /usr/local/tests/data/small_demonstration_experiment_2018_2_samples_clr_0.05.cags.json.gz \
    --temp-folder /scratch

  # Make sure the output files exist
  [[ -s test-experiment-collection.hdf5 ]]
}
