import os
import requests
import copy
import re
import pandas as pd
from ete3 import Tree
import numpy as np
import json

family_dict = None
full_tree = None

def extract_full_glottolog_tree():
    #code adapted from gerhard jaeger
    with open(raw_tree_path()) as f:
        raw = f.readlines()
    trees = []
    # each line is a tree. bring in proper format and read with ete3
    for i, ln in enumerate(raw):
        ln = ln.strip()
        ln = re.sub(r"\'[A-Z][^[]*\[", "[", ln)
        ln = re.sub(r"\][^']*\'", "]", ln)
        ln = re.sub(r"\[|\]", "", ln)
        ln = ln.replace(":1", "")
        trees.append(Tree(ln, format=1))
    # place all trees below a single root
    glot = Tree()
    for t in trees:
        glot.add_child(t)

    #insert missing, i.e. isolated languages below the root
    tTaxa = [nd.name for nd in glot.traverse() if nd.name != '']
    fn = "resources/glottolog/languages.csv"
    languages = pd.read_csv(fn)
    gTaxa = languages.Glottocode.values
    for taxon in gTaxa:
        if taxon not in tTaxa:
            glot.add_child(name=taxon)

    #if there is a inner node with a name (i.e. corresponds to a language),the name of this node is removed
    # and a child(i.e.leaf) with this name is inserted
    nonLeaves = [nd.name for nd in glot.traverse() if nd.name != '' and not nd.is_leaf()]
    for i, nm in enumerate(nonLeaves):
        nd = glot & nm
        nd.name = ''
        nd.add_child(name=nm)

    # only keep languages which are listed in languages.csv
    gTaxa = np.intersect1d(gTaxa, glot.get_leaf_names())
    glot.prune([glot&x for x in gTaxa])

    glot.write(outfile = "resources/glottolog/glottolog.tre", format=9)
    global full_tree
    full_tree = glot


def load_families():
    global family_dict
    if family_dict is not None:
        return True
    languages_path = "resources/glottolog/languages.csv"
    if not os.path.isfile(languages_path):
        return False
    languages_df = pd.read_csv(languages_path)
    family_dict = {}
    for index, row in languages_df.iterrows():
        glottocode = row["Glottocode"]
        family_id = str(row["Family_ID"])
        if family_id == "nan":
            family_dict[glottocode] = "ISOLATE"
        else:
            family_dict[glottocode] = family_id
    return True

def load_full_tree():
    global full_tree
    if full_tree is not None:
        return True
    if not os.path.isfile("resources/glottolog/glottolog.tre"):
        return False
    full_tree = Tree("resources/glottolog/glottolog.tre", format=9)
    return True

def get_family(glottocode):
    r = load_families()
    if not r:
        return None
    if glottocode not in family_dict:
        return "UNKNOWN"
    else:
        return family_dict[glottocode]


def get_families(glottocodes):
    r = load_families()
    family_ids = set()
    if not r:
        return family_ids
    for glottocode in glottocodes:
        if glottocode not in family_dict:
            family_ids.add("UNKNOWN")
            continue
        family_ids.add(family_dict[glottocode])
    return family_ids

def split_families(glottocodes):
    r = load_families()
    families = {}
    if not r:
        return families
    for glottocode in glottocodes:
        if glottocode not in family_dict:
            family_id = "UNKNOWN"
        else:
            family_id = family_dict[glottocode]
        if family_id in families:
            families[family_id].append(glottocode)
        else:
            families[family_id] = [glottocode]
    return families



def get_tree(glottocodes, languages):
    r = load_full_tree()
    if not r:
        return None
    tree = copy.deepcopy(full_tree)
    try:
        tree.prune([tree&glottocode for glottocode in glottocodes])
    except: #node not found due to wrong / deprecated glottocodes
        return None
    for leaf in tree.iter_leaves():
        leaf.add_features(new = False)
    for leaf in tree.iter_leaves():
        if leaf.new:
            continue
        ids = [languages[i] for i in range(len(glottocodes)) if glottocodes[i] == leaf.name]
        if len(ids) == 1:
            leaf.name = ids[0]
        else:
            leaf.name = ""
            for i in ids:
                leaf.add_child(name = i)
            for child in leaf.children:
                child.add_features(new = True)
    return tree
