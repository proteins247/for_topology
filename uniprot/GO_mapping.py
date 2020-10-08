#!/usr/bin/env python
"""
GO_mapping.py

Parse uniprot for GO data
"""

import re
import xml.etree.ElementTree as ET
import collections
import pandas as pd

datapath = "uniprot-proteome%3AUP000000318.xml"

tree = ET.parse(datapath)

namespace_prefix = '{http://uniprot.org/uniprot}'

# bnum aka ordered locus name (OLN)
bnum_re = re.compile("(b[0-9]+)")
go_re = re.compile("(?<=C:).*")
refseq1_re = re.compile("NP.*")
refseq2_re = re.compile("YP.*")

gather = collections.defaultdict(list)
for i, entry in enumerate(tree.iter(namespace_prefix + "entry")):
    # Get gene name, alternative names, and oln (bnum)
    accession_values = []
    for accession in entry.findall(namespace_prefix + "accession"):
        accession_values.append(accession.text)
    gather['accession'].append(accession_values)

    gene = entry.find(namespace_prefix + "gene")
    gene_name = None
    synonym_names = []
    bnum = None
    for name in gene.findall(namespace_prefix + "name"):
        if name.attrib['type'] == 'primary':
            gene_name = name.text
        elif name.attrib['type'] == 'synonym':
            synonym_names.append(name.text)
        elif name.attrib['type'] == 'ordered locus':
            if bnum_re.findall(name.text):
                bnum = name.text
    gather['name'].append(gene_name)
    gather['synonyms'].append(synonym_names)
    gather['oln'].append(bnum)

    GO_terms = []
    refseq1 = ""
    refseq2 = ""
    for reference in entry.findall(namespace_prefix + "dbReference"):
        if reference.attrib['type'] == 'GO':
            for prop in reference.findall(namespace_prefix + "property"):
                if prop.attrib['type'] == "term":
                    result = go_re.findall(prop.attrib['value'])
                    if result:
                        GO_terms.append(result[0])
        elif reference.attrib['type'] == "RefSeq":
            if refseq1_re.findall(reference.attrib['id']):
                refseq1 = reference.attrib['id']
            elif refseq2_re.findall(reference.attrib['id']):
                refseq2 = reference.attrib['id']
    gather['GO'].append(GO_terms)
    if refseq1 == "":
        gather['refseq'].append(refseq2)
    else:
        gather['refseq'].append(refseq1)

    if bnum is None:
        # Some entries will not have bnum or names, but they are not
        #   relevant?
        print("No bnum:", i, gene_name, synonym_names, refseq1)

    locations = []
    for comment in entry.findall(namespace_prefix + "comment"):
        if comment.attrib['type'] == 'subcellular location':
            subcellular_location = comment.find(namespace_prefix + "subcellularLocation")
            for location in subcellular_location.findall(namespace_prefix + "location"):
                locations.append(location.text)

    keywords = []
    for keyword in entry.findall(namespace_prefix + "keyword"):
        keywords.append(keyword.text)
    gather['keywords'].append(keywords)

dataframe = pd.DataFrame(gather)
dataframe.to_pickle("uniprot_GO+keyword.pkl")
dataframe.to_csv("uniprot_GO+keyword.csv")
print(dataframe.shape)
no_bnum = dataframe['oln'].apply(pd.isna)
dataframe = dataframe[~no_bnum]
print(dataframe.shape)
