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


def run_experiments(bn, language_set, name, num_ids, use_epitran, redo):
    if name.startswith("filter"):
        mode = "synsetfilter"
    else:
        mode = "conceptlist"


    if use_epitran:
        full_name = name + "_epitran"
    else:
        full_name = name


    bn = BabelNet.getInstance()

    assert(mode in ["conceptlist", "synsetfilter"])

    base_dir = os.path.join("results", language_set)

    if mode == "conceptlist":
        conceptlist_dir = os.path.join("resources", "conceptlist")
        conceptlist_path = os.path.join(conceptlist_dir, name + ".tsv")
        assert(os.path.isfile(conceptlist_path))

    if mode == "synsetfilter":
        synsetfilter_path = os.path.join("resources", "synsetfilter.tsv")
        ranking_dir = os.path.join(base_dir, "ranking")
        if not os.path.isdir(ranking_dir):
            os.makedirs(ranking_dir)
        ranking_path = os.path.join(ranking_dir, full_name + ".tsv")
        full_name += "_" + str(num_ids)


    babelids_dir = os.path.join(base_dir, "babelids")
    babelids_path = os.path.join(babelids_dir, full_name + "_babelids.txt")
    wordlist_dir = os.path.join(base_dir, "wordlist")
    wordlist_path = os.path.join(wordlist_dir, full_name + "_wordlist.tsv")
    wordlist_cognate_dir = os.path.join(base_dir, "wordlist_cognate")
    wordlist_cognate_path = os.path.join(wordlist_cognate_dir, full_name + "_wordlist_cognate.tsv")
    glottolog_tree_dir = os.path.join(base_dir, "glottolog_tree")
    glottolog_tree_path = os.path.join(glottolog_tree_dir, full_name + "_glottolog.tree")
    msa_dir = os.path.join(base_dir, "msa")
    bin_msa_path = os.path.join(msa_dir, full_name + ".bin.phy")
    plots_dir = os.path.join(base_dir, "plots")
    sparsity_plots_dir = os.path.join(plots_dir, "sparsity")
    sparsity_plot_path = os.path.join(plots_dir, "sparsity",  full_name + "_sparsity.png")
    statistics_dir = os.path.join(base_dir, "statistics")
    statistics_path = os.path.join(statistics_dir, full_name + "_statistics.json")
    raxml_dir = os.path.join(base_dir, "raxml", full_name)
    raxml_prefix = os.path.join(raxml_dir, "inference")
    pythia_dir = os.path.join(base_dir, "pythia")
    pythia_prefix = os.path.join(pythia_dir, full_name)
    best_tree_path = raxml_prefix + ".raxml.bestTree"

    for d in [babelids_dir, wordlist_dir, wordlist_cognate_dir, msa_dir, glottolog_tree_dir, plots_dir, sparsity_plots_dir, statistics_dir, raxml_dir, pythia_dir]:
        if not os.path.isdir(d):
            os.makedirs(d)


    langs = pipeline.get_languages(language_set) 
    if use_epitran:
        epitran_instances = util.get_epitran_instances(langs)
        #epitran_instances = [None for l in range(len(langs))]
    else:
        epitran_instances = [None for l in range(len(langs))] 
    epitran_langs = [lang for l, lang in enumerate(langs) if epitran_instances[l] is not None]

    #pipeline.languages_statistics(langs, epitran_instances)

    if mode == "conceptlist":
        pipeline.babelids_from_conceptlist(bn, conceptlist_path, babelids_path, langs, epitran_langs, redo)

    if mode == "synsetfilter":
        pipeline.filter_synsets(bn, synsetfilter_path)
        pipeline.ranking_from_synsetfilter(bn, synsetfilter_path, ranking_path, langs, epitran_langs, redo)
        pipeline.babelids_from_ranking(bn, ranking_path, babelids_path, num_ids, redo)


    if os.path.isfile(statistics_path):
        with open(statistics_path) as json_file:
            stat_dict = json.load(json_file)
    else:
        stat_dict = {}
        stat_dict["name"] = full_name
        stat_dict["num_ids"] = num_ids
        stat_dict["use_epitran"] = use_epitran

    pipeline.generate_wordlist(bn, babelids_path, wordlist_path, langs, epitran_instances, redo)
    pipeline.detect_cognates(wordlist_path, wordlist_cognate_path, redo)

    cd = CategoricalData.from_edictor_tsv(wordlist_cognate_path)
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


redo = False
bn = BabelNet.getInstance()
for language_set in ["all", "main", "ielex"]:
    for name in ["swadesh100", "swadesh200", "core-wordnet"]:
        for use_epitran in [True, False]:
            run_experiments(bn, language_set, name, float("nan"), use_epitran, redo)
    name = "filter"
    for num_ids in [100, 200, 5000]:
        for use_epitran in [True, False]:
            run_experiments(bn, language_set, name, num_ids, use_epitran, redo)

