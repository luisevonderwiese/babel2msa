import os
import re
import numpy as np
import copy
from pyglottolog import Glottolog
from ete3 import Tree


class GlottologWrapper:
    def __init__(self, glottolog_path = "../glottolog"):
        self.glottolog = Glottolog(glottolog_path)
        self.full_tree_path = os.path.join(glottolog_path, "glottolog.tre")
        if not os.path.isfile(self.full_tree_path):
            self.extract_full_glottolog_tree()
        self.full_tree = Tree(self.full_tree_path)



    def extract_full_glottolog_tree(self):
        print("Extracting Glottolog Tree ... this might take a while ...")
        #code adapted from gerhard jaeger
        raw = self.glottolog.newick_tree()
        trees = []
        # each line is a tree. bring in proper format and read with ete3
        for i, ln in enumerate(raw.split("\n")):
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
        gTaxa = [languoid.glottocode for languoid in self.glottolog.languoids(exclude_pseudo_families=True)]
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

        glot.write(outfile = self.full_tree_path, format=9)


    def get_tree(self, glottocodes, languages):
        tree = copy.deepcopy(self.full_tree)
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

