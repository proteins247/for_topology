#!/usr/bin/env python
"""
1-process_genome.py

Get genes, both dna and aa sequences

author : Victor Zhao

# Copyright (C) 2020 Victor Zhao

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

"""
import sys
import argparse
import re
import collections
import pandas as pd
from IPython import start_ipython

from Bio import SeqIO

ITEMS_RE = re.compile(r"\[.*?\]")
# GO_TO_EXCLUDE = [
#     'integral component',
#     'intrinsic component',
#     'periplasmic',
#     'anchored component of external side of plasma membrane',
#     'anchored component of cell outer membrane',
#     'external side of cell outer membrane',
#     'extracellular region'
# ]


# def check_inclusion_go(go_terms):
#     """Whether the set of GO terms (for a gene) means
#     it should be included or excluded.
#     """
#     for exclude_term in GO_TO_EXCLUDE:
#         for term in go_terms:
#             if exclude_term in term:
#                 return False
#     return True


def parse_description(description_str):
    items = ITEMS_RE.findall(description_str)
    results = {}
    for item in items:
        key, val = item.strip('[]').split('=')
        if key not in results:
            results[key] = val
        else:
            print("Duplicate key", key, val)
            print("in", description_str)
    return results


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('cds_file')
    parser.add_argument('translated_file')
    parser.add_argument('--ipython', action='store_true',
                        help="start ipython at the end of the script")
    return parser.parse_args()


def main(args):
    # Build up gene dataframe
    with open(args.cds_file) as f:
        genes = SeqIO.parse(f, "fasta")
        gather = collections.defaultdict(list)
        for gene in genes:
            annotations = parse_description(gene.description)
            if 'protein_id' not in annotations:
                if annotations['pseudo'] != "true":
                    print(gene.description)
                continue
            gather['name'].append(annotations['gene'])
            gather['bnum'].append(annotations['locus_tag'])
            gather['xref'].append(annotations['db_xref'])
            gather['description'].append(annotations['protein'])
            gather['accession'].append(annotations['protein_id'])
            gather['seq'].append(str(gene.seq))
    df_genes = pd.DataFrame(gather)
    df_genes = df_genes.set_index("accession")
    print("Genes have been read. Number:", len(df_genes))

    # Build protein dataframe to with aa sequences
    with open(args.translated_file) as f:
        proteins = SeqIO.parse(f, "fasta")
        gather = collections.defaultdict(list)
        for protein in proteins:
            annotations = parse_description(protein.description)
            if 'protein_id' not in annotations:
                continue
            # gather['name'].append(annotations['gene'])
            gather['accession'].append(annotations['protein_id'])
            gather['aaseq'].append(str(protein.seq))

    df_proteins = pd.DataFrame(gather)
    df_proteins = df_proteins.set_index("accession")
    print("CDS have been read. Number:", len(df_proteins))

    merged = pd.merge(df_genes, df_proteins, on='accession')
    print("Genes and CDS matched. Number:", len(merged))
    # I found that inner, left, or right join will yield the same result
    # 4242 genes
    # I also confirmed that all aa sequences match the length of
    # corresponding nucleotide sequences (3*aa_len + 3).

    merged.to_pickle("1_genes.pkl")
    print("Dataframe saved.")

    if args.ipython:
        sys.argv = [sys.argv[0]]
        start_ipython(user_ns=dict(locals(), **globals()))
    return


if __name__ == "__main__":
    sys.exit(main(parse_args()))
