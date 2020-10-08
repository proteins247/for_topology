#!/usr/bin/env python


import re
import glob
import collections
import pickle as pkl
import pandas as pd


enrichment_files = sorted(glob.glob("data/*rare_enrichment.dat"))

gene_name_re = re.compile(r"(\w+)_rare_enrichment.dat")

to_concat = []
for enrichment_file in enrichment_files:
    data = collections.defaultdict(list)
    gene_name = gene_name_re.findall(enrichment_file)[0]
    free_energy_file = "data/{}_free_energy.dat".format(gene_name)
    print(gene_name)
    with open(enrichment_file) as f:
        for line in f:
            if line[0] == "#":
                continue
            items = line.split()
            if len(items) != 3:
                continue
            data["name"].append(gene_name)
            data["center"].append(int(items[0]))
            data["enriched"].append(float(items[1]))
            data["p_value"].append(float(items[2]))

    dataframe = pd.DataFrame(data)
    energies = pd.read_csv(free_energy_file, sep="\s+", header=None,
                           names=["center", "F_L", "F_L_constrained"], comment='#')
    dataframe = pd.merge(dataframe, energies, on='center')
    to_concat.append(dataframe)

dataframe = pd.concat(to_concat, ignore_index=True)
dataframe.to_pickle("gene_enrichment+free_energy.pkl")
