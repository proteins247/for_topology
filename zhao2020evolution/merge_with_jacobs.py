#!/usr/bin/env python
"""
merge_with_jacobs.py

Calculate contact orders of substructures preceding rare codons.

issues:
there isn't such thing as a single offset.
insertions occur, that's ok.
should differences in aligned sequences be mapped?

accept that PDB number will be weird
https://proteopedia.org/wiki/index.php/Unusual_sequence_numbering



"""
import sys
import argparse
import collections
import gzip
import glob
import pickle as pkl
import pandas as pd
import numpy as np
np.random.seed(23)
import scipy.signal
import scipy.stats
from IPython import start_ipython

from Bio import pairwise2
import Bio.PDB
import Bio.PDB.Polypeptide
PDB_PARSER = Bio.PDB.PDBParser()
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

# Will's gene enrichment data
ENRICHMENT = "../jacobs2017evidence/gene_enrichment+free_energy.pkl"
# My database of genes and corresponding structures
STRUCTURES = "3_genes_and_downloaded_structures_best.pkl"
STRUCTURE_DIR = "structures/"
# allowed_aa_names = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS",
#                     "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL",
#                     "TRP", "TYR", "MSE"]


def expand_df(df, column):
    """Take a column of lists; give each value a row in the dataframe.

    """
    expanded2 = pd.DataFrame({
        col: np.repeat(df[col].values, df[column].str.len())
        for col in df.columns.drop(column)}
    ).assign(**{column: list(np.concatenate(df[column].values))})
    return expanded2


def calculate_long_range_order(residues, radius=7.5, exclusion=12):
    residue_resid_d = {residue: resid for resid, residue in residues}
    if residues is None:
        return None
    CAs = []
    for resid, residue in residues:
        try:
            atom = [a for a in residue.get_atoms() if a.name == 'CA'][0]
        except IndexError:
            print("Strange, residue {} lacks atom CA".format(residue.full_id))
            continue
        CAs.append(atom)
    neighbor_search = Bio.PDB.NeighborSearch(CAs)
    contact_pairs = neighbor_search.search_all(radius, level='R')

    contact_count = 0
    for contact in contact_pairs:
        resid1 = residue_resid_d[contact[0]]
        resid2 = residue_resid_d[contact[1]]
        if np.abs(resid2 - resid1) > exclusion:
            contact_count += 1
    number_residues = residues[-1][0] - residues[0][0] + 1
    return contact_count / number_residues


def get_substructure_long_range_orders(residues, locations, unextruded_length,
                                       exclusion=12):
    if residues is None:
        return None
    long_range_orders = []
    for location in locations:  # could be []
        cutoff = location - unextruded_length
        sub_residues = [
            (resid, residue) for resid, residue in residues
            if resid < cutoff]
        long_range_orders.append(calculate_long_range_order(sub_residues, exclusion=exclusion))
    return long_range_orders


def calculate_contact_order(contact_map):
    CO = 0
    for i, row in enumerate(contact_map):
        for j, contact in enumerate(row):
            if contact:         # == 1
                CO += np.abs(j - i)
    CO /= len(contact_map) * np.sum(contact_map)
    return CO


def calculate_contact_map(residues, radius=4.5, exclusion=1):
    residue_resid_d = {residue: resid for resid, residue in residues}
    atoms = []
    for resid, residue in residues:
        atoms.extend(residue.get_atoms())
    try:
        neighbor_search = Bio.PDB.NeighborSearch(atoms)
    except IndexError:
        print(residues[0][1].get_full_id()[0])
        raise
    contact_pairs = neighbor_search.search_all(radius, level='R')

    # what should the size of our contact matrix be?
    # Not the number of residues, since we should account for gaps in the structure
    # And not the max residue id, since perhaps we are missing the
    # N-terminus in some structures
    min_resid = min(residue_resid_d.values())
    map_size = max(residue_resid_d.values()) - min(residue_resid_d.values()) + 1
    contact_map = np.zeros([map_size, map_size])
    for contact in contact_pairs:
        resid1 = residue_resid_d[contact[0]]
        resid2 = residue_resid_d[contact[1]]
        if np.abs(resid2 - resid1) > exclusion:
            contact_map[resid1 - min_resid][resid2 - min_resid] = 1
    return contact_map + contact_map.T


def get_contact_order(residues, exclusion=1):
    if residues is None:
        return None
    contact_map = calculate_contact_map(residues, exclusion=exclusion)
    contact_order = calculate_contact_order(contact_map)
    return contact_order


def get_substructure_contact_orders(residues, locations, unextruded_length,
                                    exclusion=1):
    if residues is None:
        return None
    contact_orders = []
    for location in locations:  # could be []
        cutoff = location - unextruded_length
        sub_residues = [
            (resid, residue) for resid, residue in residues
            if resid < cutoff]
        contact_orders.append(get_contact_order(sub_residues, exclusion=exclusion))
    return contact_orders


def align_residues_to_sequence(residues, aaseq):
    """Return list of (resid, residue) where resid is the position of
    residue in aaseq.

    # Offset: add offset to actual PDB resid to get sequence resid
    """
    peptide = Bio.PDB.Polypeptide.Polypeptide(residues)
    peptide_seq = str(peptide.get_sequence())
    # returns a list of alignments, take top
    alignment = pairwise2.align.globalms(aaseq, peptide_seq, 2, -1, -1, -0.1,
                                         penalize_end_gaps=False)[0]
    # (sorry, i reversed aaseq and peptide_seq arguments)

    # Would be good to check that alignment is good
    score = alignment[2]
    if score / len(aaseq) < 0.85:
        print("WARNING: score / len(aaseq) < 0.85")

    # the premise is that there are no gaps for aaseq, and that
    # aaseq[0] is 1.
    aligned = []                # to be returned
    aaseq_aligned = alignment[0]
    peptide_seq_aligned = alignment[1]
    aaseq_index = 1
    peptide_seq_index = 0
    for letter_a, letter_p in zip(aaseq_aligned, peptide_seq_aligned):
        if letter_a == letter_p:
            aligned.append((aaseq_index, residues[peptide_seq_index]))
            aaseq_index += 1
            peptide_seq_index += 1
        elif letter_a == '-':
            # print("Insertion?", residues[0].full_id)
            peptide_seq_index += 1
        elif letter_p == '-':
            aaseq_index += 1
        elif letter_a != letter_p:
            # Due to mutation?
            aaseq_index += 1
            peptide_seq_index += 1
        else:
            raise Exception
    # offset = aligned[0][0] - aligned[0][1].id[1]
    # offset2 = aligned[-1][0] - aligned[-1][1].id[1]
    # if offset != offset2:
    #     print("Mismatched offsets", aligned[0][1].full_id)
    return aligned


def find_enriched_regions(name, enrichment, enrichment_threshold,
                          pvalue_threshold, excluded_length=80):
    """
    For now, just find the first peak
    """
    data = enrichment.loc[name].query("center > %d" % excluded_length)
    peaks = scipy.signal.find_peaks(data["enriched"], enrichment_threshold)[0]
    # Note that peaks are >= enrichment_threshold
    qualified_peaks = []
    for peak in peaks:
        region = slice(peak - 5, peak + 6)
        pvalues = data["p_value"].values[region]
        if np.any(pvalues < pvalue_threshold):
            qualified_peaks.append(peak)
    return data["center"].values[qualified_peaks]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ipython', action='store_true',
                        help="start ipython at the end of the script")
    parser.add_argument("--enrichment-threshold", type=float, default=0.75)
    parser.add_argument("--pvalue-threshold", type=float, default=0.01)
    parser.add_argument("--unextruded-length", type=int, default=30)
    parser.add_argument("--nterm-exclusion", type=int, default=80)
    parser.add_argument("--random-count", type=int, default=5)
    parser.add_argument("--co-exclusion", type=int, default=1)
    return parser.parse_args()


def main(args):
    # Get genes and file paths
    names = pd.read_pickle(STRUCTURES)
    # We don't need the SCOP data:
    names = names.groupby(level='name').head(1)
    names = names.drop(columns=['domain'])

    # Will's enrichment data
    enrichment = pd.read_pickle(ENRICHMENT)
    # Now, Will used different names sometimes
    # in particular, we want to change
    #   'bioD' -> 'bioD1'
    #   'umpH' -> 'nagD'
    enrichment = enrichment.set_index("name", drop=False)
    enrichment.loc['bioD', "name"] = 'bioD1'
    enrichment.loc['umpH', "name"] = 'nagD'

    # Sometimes Will used different gene names
    names['all_names'] = [[name] + synonyms for name, synonyms in
                          zip(names['name'], names['synonyms'])]
    expanded = expand_df(names, 'all_names')
    expanded = expanded.set_index(["all_names", "name"], drop=False)

    # Now for some reason, there are duplicates due to synonyms
    # I'm going to take care of the 2 that conflict
    # fabG, mgsA
    expanded = expanded.drop([('fabG', 'accC'), ('mgsA', 'rarA')])

    # Keep only those with enrichment data
    will_gene_names = [s.lower() for s in enrichment['name'].unique()]
    has_enrichment = [name.lower() in will_gene_names
                      for name in expanded.index.get_level_values('all_names')]
    expanded = expanded[has_enrichment].copy()
    # At 511 genes, exactly the number will has

    expanded["locations"] = expanded["all_names"].apply(
        find_enriched_regions,
        args=(enrichment, args.enrichment_threshold, args.pvalue_threshold),
        excluded_length=args.nterm_exclusion)
    # Note that these are locations along the gene sequence

    # Whatever further analysis:

    if args.ipython:
        sys.argv = [sys.argv[0]]
        start_ipython(user_ns=dict(locals(), **globals()))
    return 0


if __name__ == "__main__":
    sys.exit(main(parse_args()))
