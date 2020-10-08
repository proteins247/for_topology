#!/usr/bin/env python
"""
2-match_to_structures.py

Find structures for genes.

We filter at this stage too.

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
import os
import re
import sys
import argparse
import requests
import pandas as pd
import time
from IPython import start_ipython
import xml.etree.ElementTree as ET
import concurrent.futures
import numpy as np
from Bio import pairwise2


SCOP_PATH = "../dir.cla.scope.2.07-stable.txt"
RCSB_URL = "https://www.rcsb.org/pdb/rest/postBLAST"


def expand_df(df, column):
    """Take a column of lists; give each value a row in the dataframe.

    """
    expanded2 = pd.DataFrame({
        col: np.repeat(df[col].values, df[column].str.len())
        for col in df.columns.drop(column)}
    ).assign(**{column: list(np.concatenate(df[column].values))})
    return expanded2


def assess_resolution(candidates):
    """Add resolution to candidates

    Returns None if candidates is None
    Returns (None, None, None) if candidates is []
    Returns best of candidates.
    """
    if not candidates:
        return candidates       # None or []
    pdbids = ','.join([c[0] for c in candidates])
    print("Querying resolution of", pdbids)
    url = "http://www.rcsb.org/pdb/rest/describePDB?structureId={}".format(pdbids)
    status = 0
    tries = 0
    while status != 200 and tries < 5:
        try:
            response = requests.get(url)
        except:
            print("Assess resolution. Exception")
            return None
        status = response.status_code
        tries += 1
        time.sleep(1)
    if status != 200:
        print("assess_resolution failed for", candidates)
        return None
    assessed = []
    try:
        root = ET.fromstring(response.text)
    except:
        return None
    for pdb, candidate in zip(root.getchildren(), candidates):
        assert(pdb.attrib["structureId"].upper() == candidate[0].upper())
        if pdb.attrib["expMethod"].upper() == "SOLUTION NMR":
            # Give NMR a resolution of 2
            resolution = 2.0
        else:
            try:
                resolution = float(pdb.attrib["resolution"])
            except KeyError:
                resolution = 100
        assessed.append(
            candidate + (resolution,))
    return assessed


def candidates_annotated_by_scop(candidates, scop_data):
    """Mark if candidates have been annotated by SCOP.
    """
    if not candidates:
        return candidates       # either None or []
    modified_candidates = []
    for candidate in candidates:
        pdbid = candidate[0]
        modified_candidates.append(
            candidate + (pdbid.lower() in scop_data['pdbid'].values,))
    return modified_candidates


def parse_response(response, original_sequence, identity_threshold):
    """For a RCSB response, return candidates with identities of at
    least length * identity_threshold.

    Returns None if response is XML without results due to 503
    Returns [] if there aren't any hits
    Returns a list of (pdbid, chain(s), # identical residues)

    We count identity residues by doing alignment. Why? Because the
    blast search masks low complexity regions, and I can't get rid of
    that.
    """
    if not response:
        return None
    try:
        root = ET.fromstring(response.text)
    except:
        return None
    try:
        hits = root.find("BlastOutput_iterations").find("Iteration").find("Iteration_hits")
    except AttributeError:
        print("AttributeError. Failed to get Iteration_hits")
        print(response.text)
        return None
    if hits is None:
        return []
    query_length = int(root.find("BlastOutput_query-len").text)
    hit_threshold = query_length * identity_threshold
    to_assess = []
    below_count = 0
    for hit in hits.getchildren():
        definition = hit.find("Hit_def").text
        hit_data = hit.find("Hit_hsps").find("Hsp")
        hit_seq = hit_data.find("Hsp_hseq").text
        identity_residues = pairwise2.align.globalxx(original_sequence, hit_seq, score_only=True)
        # identity_residues = int(hit_data.find("Hsp_identity").text)
        if identity_residues < hit_threshold:
            below_count += 1
            if below_count >= 3:
                break
            continue
        else:
            below_count = 0
        pdbid, _, chains = definition.split('|')[0].split(':')
        to_assess.append((pdbid, chains, identity_residues))
    return to_assess


def query_rcsb(sequence):
    """For a particular sequence, request matching structures.

    Sometimes this fails (too many queries?)
    """
    print('Querying "{}..."'.format(sequence[:4]))
    data = {
        "sequence": sequence,
        "eCutOff": 10.0,
        # "maskLowComplexity": "no",
        "matrix": "BLOSUM62",
        "outputFormat": "XML"
    }
    status = 0
    tries = 0
    while status != 200 and tries < 5:
        try:
            ret = requests.post(RCSB_URL, data=data)
        except:
            print("Query rcsb. Exception")
            return None
        status = ret.status_code
        tries += 1
        time.sleep(1)
    if status != 200:
        print('Failed querying "{}"'.format(sequence[:10]))
        return None
    return ret


def parse_scop():
    scop_data = pd.read_csv(
        "../SCOP/dir.cla.scope.2.07-stable.txt", comment='#', header=None, sep=r'\s+',
        names=["sid", "pdbid", "domain", "sccs", "px", "sunids"])
    return scop_data.set_index("pdbid", drop=False)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('uniprot_GO')
    parser.add_argument('--ipython', action='store_true',
                        help="start ipython at the end of the script")
    parser.add_argument('--skip-reevaluation', action='store_true')
    parser.add_argument('--identity-threshold', default=0.85, type=float)
    return parser.parse_args()


def main(args):
    # SCOP data is needed for checking if annotated
    scop_data = parse_scop()

    # Perform queries if no existing file
    if not os.path.isfile("2_genes_with_structures.pkl"):
        data = pd.read_pickle("1_genes.pkl")

        # Threads speeds this up, but maybe causes some queries to fail?
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
            responses = executor.map(query_rcsb, data['aaseq'])

        all_candidates = [parse_response(response, aaseq, args.identity_threshold)
                          for response, aaseq in zip(responses, data['aaseq'])]
        all_candidates = [candidates_annotated_by_scop(candidates, scop_data)
                          for candidates in all_candidates]

        # structure = [assess_resolution(candidates)
        #              for candidates in all_candidates]
        # print(structure)
        # sys.exit(0)
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
            assessed = executor.map(assess_resolution, all_candidates)
        data["structure"] = list(assessed)
    else:
        data = pd.read_pickle("2_genes_with_structures.pkl")
        print("Opened existing file.")
        if not args.skip_reevaluation:
            # Because sometimes queries seem to fail, we'll check. Look
            # at structures with None
            for idx, sequence, structure in zip(data.index, data['aaseq'], data['structure']):
                if structure is None:
                    response = query_rcsb(sequence)
                    candidates = parse_response(response, sequence, args.identity_threshold)
                    candidates = candidates_annotated_by_scop(candidates, scop_data)
                    assessed = assess_resolution(candidates)
                    data.at[idx, "structure"] = assessed
                    print("Updated", idx)
    print("Writing file.")
    data.to_pickle("2_genes_with_structures.pkl")

    # Add in information. GO merge
    # ----------------------------------------

    data = data[data['structure'].apply(len) != 0].copy()
    print("Removed genes without structures. New size:", len(data))

    data["aalen"] = data['aaseq'].apply(len)

    # Merging on ncbi accession without version nujmbers works best
    data = data.reset_index()
    data['accession_s'] = data['accession'].apply(lambda x: x.split('.')[0])
    uniprot_data = pd.read_pickle(args.uniprot_GO)
    uniprot_data['refseq_s'] = uniprot_data['refseq'].apply(lambda x: x.split('.')[0])
    merged = pd.merge(data.reset_index(), uniprot_data,
                      how='left', left_on="accession_s", right_on="refseq_s")
    merged.drop(inplace=True, columns=["xref", "accession_x", "accession_y", 'index',
                                       'accession_s', 'refseq_s', "name_x", 'oln'])
    merged['name'] = merged['name_y']
    merged.drop(inplace=True, columns=['name_y'])
    data = merged
    print("Genes+CDS matched to uniprot. Number:", len(merged))
    print(data.head(10))

    # SCOP merge
    # ----------------------------------------

    # Sort structures by SCOP annotated, sequence coverage, and
    # resolution. Keep top 2.
    data['structure'] = data['structure'].apply(
        sorted, key=lambda x: (~x[3], -x[2], x[4]))

    def up_to_top2(structures):
        to_keep = structures[:1]
        if len(structures) < 2:
            return to_keep
        else:
            if structures[0][3] == structures[1][3]:
                to_keep.append(structures[1])
            return to_keep

    data['structure'] = data['structure'].apply(up_to_top2)
    data = expand_df(data, 'structure')
    # print("Number with SCOP annotations:",
    #       np.sum(data['structure'].apply(lambda x: x[0][3])))
    # print("  Note, non-annotated are being kept too.")

    # Merge in SCOP data
    data['pdbid'] = data['structure'].apply(lambda x: x[0].lower())
    data['chain'] = data['structure'].apply(lambda x: x[1])
    scop_to_merge = scop_data.reset_index(drop=True)[["pdbid", "domain"]]
    merged_scop = pd.merge(
        data, scop_to_merge, how='left', left_on="pdbid", right_on="pdbid")

    def check_row_relevant(subdf):
        """Only return relevant rows.
        """
        chains = subdf.iloc[0]['structure'][1].split(',')  # as identified by RCSB
        scop_chains = []                                   # what SCOP annotated
        for value in subdf['domain']:
            if pd.isna(value):
                continue
            scop_chains.append(value.split(':')[0])
        if not set(chains).intersection(scop_chains):
            # no overlap because SCOP annotation is for chain not
            #   relevant to our gene.
            ret = subdf.head(1).copy()
            ret["domain"] = np.nan
            return ret
        elif not scop_chains:
            # Because there is no domain annotation
            truth_array = [True] * len(subdf)
            if len(subdf) > 1:
                raise Exception("unexpected")
        else:
            truth_array = [c in chains for c in scop_chains]
        return subdf[truth_array]

    merged_scop = merged_scop.groupby(["name", "pdbid", "chain"]).apply(check_row_relevant)
    # merged_scop.drop(inplace=True, columns=['chain'])

    print("After merging in SCOP, # genes:", len(merged_scop['name'].unique()))
    print("After merging in SCOP, # rows:", len(merged_scop))
    # print("After merging in SCOP, # genes:", len(merged_scop['name'].unique()))
    # print("After merging in SCOP, # domains:", len(merged_scop))
    merged_scop.to_pickle("2_genes_with_structures_filtered.pkl")
    print("Saved to file")

    if args.ipython:
        sys.argv = [sys.argv[0]]
        start_ipython(user_ns=dict(locals(), **globals()))
    return


if __name__ == "__main__":
    sys.exit(main(parse_args()))
