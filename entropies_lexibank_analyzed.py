import os
import json
import traceback

from cognate import CognateData


import jpype
import jpype.imports
from jpype.types import *

# Launch the JVM
jpype.startJVM(classpath=['BabelNet-API-5.3/lib/*', 'BabelNet-API-5.3/babelnet-api-5.3.jar', 'config'])


from it.uniroma1.lcl.babelnet import BabelNet
from it.uniroma1.lcl.jlt.util import Language

import util
import pipeline
import matplotlib.pyplot as plt
language_set = "lexibank-analyzed"
domain = "families"
redo = False

base_dir = os.path.join("results", language_set + "_" + domain)
wordlist_dir = os.path.join("resources", "lexibank-analyzed_wordlists", domain)
wordlist_cognate_dir = os.path.join(base_dir, "wordlist_cognate")
plots_dir = os.path.join(base_dir, "plots")
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)
families = set()
for file_name in os.listdir(wordlist_dir):
    families.add(file_name.split("_")[0])

all_entropies = []
for family in families:
    full_name = family 
    wordlist_path = os.path.join(wordlist_dir, full_name + "_wordlist.tsv")
    if not os.path.isfile(wordlist_path):
        continue
    print(full_name)
    wordlist_cognate_path = os.path.join(wordlist_cognate_dir, full_name + "_wordlist_cognate.tsv")
    try:
        pipeline.detect_cognates(wordlist_path, wordlist_cognate_path, redo)
    except Exception as e:
        traceback.print_exc()
        print(e)
        continue
    cd = CognateData.from_edictor_tsv(wordlist_cognate_path)
    if cd.num_languages() < 4:
        print(family, "too small")
        continue
    all_entropies.append(cd.bin_entropy())

plt.hist(all_entropies, bins = 20)
plt.xlabel("entropy")
plt.ylabel("num datasets")
plt.savefig(os.path.join(plots_dir, "hist_entropies.png"))
plt.clf()
plt.close()
