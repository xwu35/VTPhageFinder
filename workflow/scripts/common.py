#!/usr/bin/env python

import os
import glob
import pandas as pd

# function to get variables and sequences path
def parse_samples_and_sequences(metadata, raw_reads_dir=""):
    R1_map = {}
    R2_map = {}
    df = pd.read_table(metadata).set_index("sample", drop=False)
    samples = df.index.tolist()
    for sample in samples:
      sample_R1 = df.loc[sample, "R1"]
      sample_R2 = df.loc[sample, "R2"]
      R1_map[sample] = os.path.join(raw_reads_dir, sample_R1)
      R2_map[sample] = os.path.join(raw_reads_dir, sample_R2)
    return samples, R1_map, R2_map

