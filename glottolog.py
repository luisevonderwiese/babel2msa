import os
from pyglottolog import Glottolog
from ete3 import Tree

class GlottologWrapper:
    def __init__(self, glottolog_path = "../glottolog"):
        print("Initializing glottolog...")
        self.glottolog = Glottolog(glottolog_path)
        self.languoid_dict = self.glottolog.languoids_by_code()
        self.lgs = {lg.id: lg for lg in self.glottolog.languoids()}
        print("done")

    def get_iso(self, glottocode):
        l = self.glottolog.languoid(glottocode)
        if l:
            return l.iso
        return None

    def get_glotto(self, iso):
        if not iso in self.languoid_dict:
            return None
        return self.languoid_dict[iso].glottocode
    
    def get_glottocodes(self, iso_codes):
        all_glottocodes = []
        for code in iso_codes:
            glottocode = self.get_glotto(code)
            assert(glottocode)
            if glottocode:
                all_glottocodes.append(glottocode)
            else:
                all_glottocodes.append("")
        return all_glottocodes


    def get_tree(self, glottocodes, languages):
        family_gcs = set() 
        for gc in glottocodes:
            if self.lgs[gc].lineage:
                family_gcs.add(self.lgs[gc].lineage[0][1])
            else:
                family_gcs.add(gc)
        tree = Tree()
        for family_gc in family_gcs:
            tree_str = "(" + str(self.lgs[family_gc].newick_node(template='{l.id}').newick) + ");"
            family_tree = Tree(tree_str, format = 1)
            tree.add_child(family_tree)
        try:
            tree.prune([tree&glottocode for glottocode in glottocodes])
        except Exception as e: #node not found due to wrong / deprecated glottocodes
            print(e)
            return None

        #case that glottocodes are assigned to inner nodes
        for node in tree.traverse():
            if not node.is_leaf() and node.name in glottocodes:
                gc =  node.name
                node.name = ""
                node.add_child(name = gc)
                node.resolve_polytomy(recursive=False)

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
