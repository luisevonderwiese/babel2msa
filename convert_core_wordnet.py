import os

import jpype
import jpype.imports
from jpype.types import *

# Launch the JVM
jpype.startJVM(classpath=['BabelNet-API-5.3/lib/*', 'BabelNet-API-5.3/babelnet-api-5.3.jar', 'config'])
import util

with open(os.path.join("resources", "conceptlist", "core-wordnet_raw.txt"), "r") as infile:
    lines = infile.readlines()

concepticon_map = util.get_concepticon_map()
to_proper_pos = {"n": "NOUN", "a": "ADJ", "v": "VERB"}
with open(os.path.join("resources", "conceptlist", "core-wordnet.tsv"), "w+") as outfile:
    outfile.write("\t".join(["ID", "CONCEPTICON_ID","CONCEPTICON_GLOSS","NUMBER","ENGLISH","POS"]) + "\n")
    for ID, line in enumerate(lines):
        parts = line[:-1].split(" ")
        concept = parts[2].strip("[]")
        pos = to_proper_pos[parts[0]]
        if concept in concepticon_map:
            concepticon_id, concepticon_gloss = concepticon_map[concept]
        else:
            concepticon_id, concepticon_gloss = ("", "")
        outfile.write("\t".join(["core-wordnet-" + str(ID), str(concepticon_id), concepticon_gloss, str(ID), concept, pos]) + "\n")



