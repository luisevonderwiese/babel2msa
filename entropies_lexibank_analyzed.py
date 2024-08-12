import os
import json
import traceback

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
language_set = "lexibank-analyzed"
domain = "families"
redo = False

base_dir = os.path.join("results", language_set + "_" + domain)
wordlist_dir = os.path.join("resources", "lexibank-analyzed_wordlists", domain)
wordlist_cognate_dir = os.path.join(base_dir, "wordlist_cognate")
plots_dir = os.path.join(base_dir, "plots")

families = set()
for file_name in os.listdir(wordlist_dir):
    families.add(file_name.split("_")[0])

all_entropies = []
for family in families:
    full_name = family + "_all" 
    wordlist_path = os.path.join(wordlist_dir, full_name + "_wordlist.tsv")
    if not os.path.isfile(wordlist_path):
        continue
    print(full_name)
    wordlist_cognate_path = os.path.join(wordlist_cognate_dir, full_name + "_wordlist_cognate.tsv")


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
