import os
import json
from categorical import CategoricalData


import jpype
import jpype.imports
from jpype.types import *

# Launch the JVM
jpype.startJVM(classpath=['BabelNet-API-5.3/lib/*', 'BabelNet-API-5.3/babelnet-api-5.3.jar', 'config'])


from it.uniroma1.lcl.babelnet import BabelNet
from it.uniroma1.lcl.jlt.util import Language

import util
import pipeline
import patch_lingpy

language_set = "lexibank-analyzed"
domain = "languagelists"
redo = False

base_dir = os.path.join("results", language_set + "_" + domain)
wordlist_dir = os.path.join("resources", "lexibank-analyzed_wordlists", domain)
wordlist_cognate_dir = os.path.join(base_dir, "wordlist_cognate")
glottolog_tree_dir = os.path.join(base_dir, "glottolog_tree")
msa_dir = os.path.join(base_dir, "msa")
plots_dir = os.path.join(base_dir, "plots")
sparsity_plots_dir = os.path.join(plots_dir, "sparsity")
statistics_dir = os.path.join(base_dir, "statistics")
pythia_dir = os.path.join(base_dir, "pythia")

for d in [wordlist_cognate_dir, msa_dir, glottolog_tree_dir, plots_dir, sparsity_plots_dir, statistics_dir, pythia_dir]:
    if not os.path.isdir(d):
        os.makedirs(d)
families = set()
for file_name in os.listdir(wordlist_dir):
    families.add(file_name.split("_")[0])


for family in families:
    for conceptlist in ["swadesh100", "swadesh200", "all"]:
        full_name = family + "_" + conceptlist 
        wordlist_path = os.path.join(wordlist_dir, full_name + "_wordlist.tsv")
        if not os.path.isfile(wordlist_path):
            continue
        print(full_name)
        wordlist_cognate_path = os.path.join(wordlist_cognate_dir, full_name + "_wordlist_cognate.tsv")
        glottolog_tree_path = os.path.join(glottolog_tree_dir, full_name + "_glottolog.tree")
        bin_msa_path = os.path.join(msa_dir, full_name + ".bin.phy")
        sparsity_plot_path = os.path.join(plots_dir, "sparsity",  full_name + "_sparsity.png")
        statistics_path = os.path.join(statistics_dir, full_name + "_statistics.json")
        raxml_dir = os.path.join(base_dir, "raxml", full_name)
        raxml_prefix = os.path.join(raxml_dir, "inference")
        pythia_prefix = os.path.join(pythia_dir, full_name)
        best_tree_path = raxml_prefix + ".raxml.bestTree"

        if not os.path.isdir(raxml_dir):
            os.makedirs(raxml_dir)

        if os.path.isfile(statistics_path):
            with open(statistics_path) as json_file:
                stat_dict = json.load(json_file)
        else:
            stat_dict = {}
            stat_dict["name"] = full_name
            stat_dict["family"] = family
            stat_dict["conceptlist"] = conceptlist

        try:
            pipeline.detect_cognates(wordlist_path, wordlist_cognate_path, redo)
        except Exception as e: 
            print(e)
            continue
        cd = CategoricalData.from_edictor_tsv(wordlist_cognate_path)
        if cd.num_taxa() < 4:
            continue
        if not os.path.isfile(bin_msa_path) or redo:
            print("Writing MSA")
            cd.write_msa(bin_msa_path, "bin")
        if not os.path.isfile(glottolog_tree_path) or redo:
            print("Writing tree")
            tree = cd.get_glottolog_tree()
            tree.write(format = 1, outfile = glottolog_tree_path)

        pipeline.sparsity_plot(cd, sparsity_plot_path, redo)
        stat_dict = pipeline.statistics(wordlist_path, plots_dir, stat_dict, redo)
        stat_dict = pipeline.cognate_statistics(wordlist_cognate_path, plots_dir, stat_dict, redo)
        pipeline.run_raxmlng(bin_msa_path, "BIN+G", raxml_prefix, redo)
        if "gq_dist" not in stat_dict or redo:
            stat_dict["gq_dist"]  = util.gq_distance(glottolog_tree_path, best_tree_path)
        print("GQ distance to glottolog:", str(stat_dict["gq_dist"]))

        pipeline.run_pythia(bin_msa_path, pythia_prefix, redo)
        if "difficulty" not in stat_dict or redo:
            stat_dict["difficulty"] = util.get_difficulty(pythia_prefix)

        with open(statistics_path, "w+") as outfile:
            json.dump(stat_dict, outfile)



