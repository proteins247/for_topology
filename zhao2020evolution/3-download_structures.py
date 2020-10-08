#!/usr/bin/env python
"""
3-download_structures.py

Donwload the structures, identifying best chains.

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
import sys
import gzip
import glob
import argparse
import requests
import functools
import collections
import numpy as np
import pandas as pd
from IPython import start_ipython

import prody
from Bio import pairwise2

DOWNLOADDIR = "structures/"
PREURL = "https://files.rcsb.org/download/{}"


class memoized():
    '''Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
    (not reevaluated).
    '''
    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        if not isinstance(args, collections.Hashable):
            # uncacheable. a list, for instance.
            # better to not cache than blow up.
            return self.func(*args)
        if args in self.cache:
            return self.cache[args]
        else:
            value = self.func(*args)
            self.cache[args] = value
            return value

    def __repr__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__

    def __get__(self, obj, objtype):
        '''Support instance methods.'''
        return functools.partial(self.__call__, obj)


@memoized
def identify_best_chain(pdbid, aaseq, chain_str):
    """Determine PDB file and chain that best matches the aa sequence
    """
    files = glob.glob(DOWNLOADDIR + pdbid + "*pdb*")
    chains = chain_str.split(',')

    chain_scores = []
    for file_path in files:
        if file_path.split('.')[-1] == "gz":
            f = gzip.open(file_path, 'rt')
            filter_chains = True
        else:
            f = open(file_path, "r")
            filter_chains = False
        try:
            mol, header = prody.parsePDBStream(f, header=True)
        except ValueError as e:
            print("PDB Error: %s: %s" % (pdbid, e))
            return None
        f.close()

        for chain in mol.iterChains():
            if filter_chains and chain.getChid() not in chains:
                continue
            chain_sequence = chain.getSequence()
            # score = pairwise2.align.globalms(aaseq, chain_sequence, 2, -1, -1, -.5, score_only=True)
            # ms: https://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
            # 2 for identity, -1 for non-identity, -1 for gap, -0.5 for extending gap
            score = pairwise2.align.globalxx(aaseq, chain_sequence, score_only=True)
            # xx: no parameters
            if isinstance(score, list):
                # case [], doesn't align at all
                continue
            chain_scores.append((chain, score, file_path))
    best = sorted(chain_scores, key=lambda x: x[1], reverse=True)[0]
    best_chain = best[0]
    best_file = best[2].lstrip(DOWNLOADDIR)
    best_chain_id = best_chain.getChid()
    best_score = best[1] / float(len(aaseq))
    return best_file, best_chain_id, best_score


def download_structure(pdbid):
    filename = pdbid + ".pdb.gz"
    if os.path.isfile(DOWNLOADDIR + filename):
        print("Already downloaded", pdbid)
        return True
    if glob.glob(DOWNLOADDIR + pdbid + "*.pdb"):
        print("Already downloaded", pdbid)
        return True
    print("Downloading", pdbid)
    url = PREURL.format(filename)
    response = requests.get(url)
    if response.reason == 'OK':
        with open(DOWNLOADDIR + filename, "wb") as f:
            f.write(response.content)
        print("Saved", filename)
        return True
    else:
        print(" Failed to download")
        return False


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ipython', action='store_true',
                        help="start ipython at the end of the script")
    return parser.parse_args()


def main(args):
    failed_pdbids = []
    if not os.path.isfile("3_genes_and_downloaded_structures.pkl"):
        data = pd.read_pickle("2_genes_with_structures_filtered.pkl")

        pdbids = data['pdbid']

        filenames = []
        best_chains = []             # chain and residues
        sequence_coverages = []
        for pdbid, aaseq, chain in zip(pdbids, data['aaseq'], data['chain']):
            success = download_structure(pdbid)
            if success:
                filename, best_chain, coverage = identify_best_chain(pdbid, aaseq, chain)
                filenames.append(filename)
                best_chains.append(best_chain)
                sequence_coverages.append(coverage)
            else:
                failed_pdbids.append(pdbid)
                filenames.append(None)
                best_chains.append(None)
                sequence_coverages.append(None)

        data["filename"] = filenames
        data["best_chain"] = best_chains
        data["sequence_coverage"] = sequence_coverages
    else:
        data = pd.read_pickle("3_genes_and_downloaded_structures.pkl")
        for idx, pdbid, aaseq, chain, filename in zip(
                data.index, data['pdbid'], data['aaseq'], data['chain'], data['filename']):
            if filename is None:
                success = download_structure(pdbid)
                if success:
                    filename, best_chain, coverage = identify_best_chain(pdbid, aaseq, chain)
                    data.at[idx, "filename"] = filename
                    data.at[idx, "best_chain"] = best_chain
                    data.at[idx, "sequence_coverage"] = coverage
                else:
                    failed_pdbids.append(pdbid)
    data.to_pickle("3_genes_and_downloaded_structures.pkl")
    print("Saved genes and downloaded structures. # rows:", len(data))

    if failed_pdbids:
        print("Failed to download these pdbids:", list(set(failed_pdbids)))
        print("The failure is likely due to these structures being too large to fit\n"
              "in single PDB files. Download manually.\n"
              "https://www.rcsb.org/pages/download_features")
        print(','.join(list(set(failed_pdbids))))

    data = data[data['filename'].apply(lambda x: x is not None)]

    # ----------------------------------------
    # First, get the best structure
    # ----------------------------------------
    def best_structure(subdf, resolution_limit=4):
        """
        case 1: just a single structure
        case 2: coverages differ
        """
        first_rows = subdf.groupby(level='pdbid').head(1)
        pdbids = first_rows['pdbid'].values
        chains = first_rows['chain'].values
        if len(pdbids) == 1:
            # Only 1 pdb structure for this gene
            return subdf.xs([pdbids[0], chains[0]], level=['pdbid', 'chain'])
        assert len(pdbids) == 2
        resolutions = np.array([float(v[-1]) for v in first_rows['structure']])
        if np.any(resolutions > resolution_limit) and not np.all(resolutions > resolution_limit):
            # 1 structure is over the resolution limit, so
            # automatically choose the other
            winner_p = pdbids[np.argmin(resolutions)]
            winner_c = chains[np.argmin(resolutions)]
            return subdf.xs([winner_p, winner_c], level=['pdbid', 'chain'])
        else:
            # Choose pdbid with best coverage
            coverages = first_rows['sequence_coverage'].values
            winner_p = pdbids[np.argmax(coverages)]
            winner_c = chains[np.argmax(coverages)]
            return subdf.xs([winner_p, winner_c], level=['pdbid', 'chain'])

    data = data.groupby(level='name').apply(best_structure)
    data = data.reset_index(drop=True).set_index("name", drop=False)
    data.to_pickle("3_genes_and_downloaded_structures_best.pkl")

    print("Identified best PDB structure (out of max 2 per gene)")
    print("Saved genes and downloaded structures to '3_genes_and_downloaded_structures_best.pkl'")
    print("# rows:", len(data))

    # Also, we've identified the best chain, so we don't need to keep all SCOP entries
    def check_row_relevant(subdf):
        """For a particular gene, we've identified the best chain,
        so we want to get rid of all irrelevant SCOP annotations.

        i.e. if the best chain is A, then we don't need domain
        information for B, C, ...

        Either domain is NaN or has SCOP data. NaN means no SCOP.
        """
        assert len(subdf['best_chain'].unique()) == 1
        best_chain = subdf['best_chain'].unique()[0]
        scop_chains = []
        for value in subdf['domain']:
            if pd.isna(value):
                scop_chains.append(best_chain)
            else:
                scop_chains.append(value.split(':')[0])
        truth_array = np.array(scop_chains) == best_chain
        return subdf[truth_array]

    data = data.groupby(level='name', group_keys=False).apply(check_row_relevant)
    print("We've identified best chain, so exclude SCOP entries that don't match")
    print("Remaining # domains", len(data))
    print("Remaining # pdbfiles", len(data['pdbid'].unique()))
    print("saving")

    data.to_pickle("3_genes_and_downloaded_structures_best_filtered.pkl")

    if args.ipython:
        sys.argv = [sys.argv[0]]
        start_ipython(user_ns=dict(locals(), **globals()))
    return


if __name__ == "__main__":
    sys.exit(main(parse_args()))
