import os
from categorical import CategoricalData

import jpype
import jpype.imports
from jpype.types import *

# Launch the JVM
jpype.startJVM(classpath=['BabelNet-API-5.3/lib/*', 'BabelNet-API-5.3/babelnet-api-5.3.jar', 'config'])

import pipeline

plots_dir = os.path.join("results", "completeness_plots")
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)

for pre in ["main", "iecor"]:
    lexibank_path = os.path.join("results", "lexibank-analyzed_languagelists", "wordlist_cognate",  pre + "_swadesh100_wordlist_cognate.tsv")
    babelnet_path = os.path.join("results", pre, "wordlist_cognate", "swadesh100_epitran_wordlist_cognate.tsv")
    lexibank_cd = CategoricalData.from_edictor_tsv(lexibank_path)
    babelnet_cd = CategoricalData.from_edictor_tsv(babelnet_path)
    lexibank_glottocodes = set(lexibank_cd.glottocodes)
    babelnet_glottocodes = set(babelnet_cd.glottocodes)
    common_glottocodes = lexibank_glottocodes.intersection(babelnet_glottocodes)
    lexibank_subcd = lexibank_cd.subset_with_glottocodes(common_glottocodes)
    babelnet_subcd = babelnet_cd.subset_with_glottocodes(common_glottocodes)
    pipeline.sparsity_plot(lexibank_subcd, os.path.join(plots_dir, pre + "_lexibank_sparsity.png"), True)
    pipeline.sparsity_plot(babelnet_subcd, os.path.join(plots_dir, pre + "_babelnet_sparsity.png"), True)
