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

import matplotlib.pyplot as plt

language_set = "lexibench"
redo = False

base_dir = os.path.join("results", language_set)
wordlist_cognate_dir = os.path.join("resources", "lexibench_wordlists")
plots_dir = os.path.join(base_dir, "plots")

datasets = set()
for file_name in os.listdir(wordlist_cognate_dir):
    datasets.add(file_name.split("_")[0])

all_entropies = []
for dataset in datasets:
    full_name = dataset + "_all"
    wordlist_cognate_path = os.path.join(wordlist_cognate_dir, full_name + "_wordlist_cognate.tsv")
    if not os.path.isfile(wordlist_cognate_path):
        continue

    print(full_name)
    try:
        cd = CategoricalData.from_edictor_tsv(wordlist_cognate_path)
    except:
        continue
    if cd.num_taxa() < 4:
        continue
    all_entropies.append(cd.bin_entropy())

plt.hist(all_entropies, bins = 20)
plt.xlabel("entropy")
plt.ylabel("num datasets")
plt.savefig(os.path.join(plots_dir, "hist_entropies.png"))
plt.clf()
plt.close()



