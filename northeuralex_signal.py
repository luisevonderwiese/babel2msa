import os
import pandas as pd

import epitran
from epitran.backoff import Backoff

from ipatok import tokenise

from categorical import CategoricalData


import jpype
import jpype.imports
from jpype.types import *

# Launch the JVM
jpype.startJVM(classpath=['BabelNet-API-5.3/lib/*', 'BabelNet-API-5.3/babelnet-api-5.3.jar', 'config'])


import pipeline
import util

def get_epitran_dict():
    epitran_dict = {}
    with open(os.path.join("resources", "epitran_codes.txt"), "r") as codes_file:
        epitran_codes = codes_file.readlines()
        epitran_codes = [code[:-1] for code in epitran_codes]
    for code in epitran_codes:
        lang = code.split("-")[0]
        if not lang in epitran_dict:
            epitran_dict[lang] = []
        epitran_dict[lang].append(code)
    return epitran_dict


def get_epitran(code, epitran_dict):
    if code in epitran_dict:
        epitran_codes = epitran_dict[code]
        if len(epitran_codes) == 1:
            try:
                epi = epitran.Epitran(epitran_codes[0])
                return epi
            except Exception as e:
                pritn("epi failed")
                print(e)
                return None
        else:
            try:
                if code == "cmn":
                    epi = Backoff(epitran_codes, cedict_file=os.path.join("resources", "cedict_1_0_ts_utf-8_mdbg.txt"))
                else:
                    epi = Backoff(epitran_codes)
                return epi
            except Exception as e:
                print("backoff failed")
                print(e)
                return None
    else:
        return None



def generate_wordlists(wordlist_paths):
    df = pd.read_csv(os.path.join("resources", "northeuralex-0.9-forms.tsv"), sep = "\t")
    df = df.astype("str")
    gi_map = util.get_glotto_iso_map()
    epitran_dict = get_epitran_dict()
    epitran_instances = {}
    for glottocode in set(df["Glottocode"]):
        if glottocode in gi_map:
            print(glottocode)
            iso_code = gi_map[glottocode]
            epi = get_epitran(iso_code, epitran_dict)
            if epi is not None:
                epitran_instances[iso_code] = epi

    for path in wordlist_paths: 
        with open(path, "w+") as wordlist_file:
            wordlist_file.write("\t".join(["ID","DOCULECT","GLOTTOCODE", "ISO_CODE", "CONCEPT","CONCEPTICON_ID", "CONCEPTICON_GLOSS", "FORM", "IPA", "TOKENS"]) + "\n")

    ID = 0
    for c, row in df.iterrows():
        doculect = row["Language_ID"]
        glottocode = row["Glottocode"]
        if glottocode not in gi_map:
            continue
        iso_code = gi_map[glottocode]
        if iso_code not in epitran_instances:
            continue
        concept = row["Concept_ID"]
        print(concept)
        concepticon_id = ""
        concepticon_gloss = ""
        form = row["Word_Form"]
        ipa = row["rawIPA"]
        ipa, valid = util.process_ipa(ipa)
        if not valid:
            continue
        tokens = row["IPA"]
        with open(wordlist_paths[0], "a") as wordlist_file:
            wordlist_file.write("\t".join([str(ID), doculect, glottocode, iso_code, concept, concepticon_id, concepticon_gloss, form, ipa, tokens]) + "\n")
        tokens = tokenise(ipa, strict=False, replace=True, diphthongs=False, tones=False, unknown=False)
        tokens = list(filter(lambda token: token not in ["ʼ", "`", "´", "", " "], tokens))
        tokens = " ".join(tokens) 
        if len(tokens) > 0:
            with open(wordlist_paths[1], "a") as wordlist_file:
                wordlist_file.write("\t".join([str(ID), doculect, glottocode, iso_code, concept, concepticon_id, concepticon_gloss, form, ipa, tokens]) + "\n")
        epi = epitran_instances[iso_code]
        ipa = epi.transliterate(form)
        ipa, valid = util.process_ipa(ipa)
        if not valid:
            continue
        tokens = " ".join(tokenise(ipa))
        tokens = tokenise(ipa, strict=False, replace=True, diphthongs=False, tones=False, unknown=False)
        tokens = list(filter(lambda token: token not in ["ʼ", "`", "´", "", " "], tokens))
        tokens = " ".join(tokens)
        if len(tokens) > 0:
            with open(wordlist_paths[2], "a") as wordlist_file:
                wordlist_file.write("\t".join([str(ID), doculect, glottocode, iso_code, concept, concepticon_id, concepticon_gloss, form, ipa, tokens]) + "\n")
        ID += 1

results_dir = os.path.join("results", "northeuralex")
names = ["original", "ipatok", "epitran"]
redo = False

wordlist_dir = os.path.join(results_dir,"wordlist")
wordlist_cognate_dir = os.path.join(results_dir, "wordlist_cognate")
glottolog_tree_dir = os.path.join(results_dir, "glottolog_tree")
msa_dir = os.path.join(results_dir, "msa")
pythia_dir = os.path.join(results_dir, "pythia")

for d in [wordlist_dir, wordlist_cognate_dir, msa_dir, glottolog_tree_dir, pythia_dir]:
    if not os.path.isdir(d):
        os.makedirs(d)


wordlist_paths = [os.path.join(wordlist_dir,  name + "_wordlist.tsv") for name in names]
#generate_wordlists(wordlist_paths)


for name in names:
    wordlist_path = os.path.join(wordlist_dir,  name + "_wordlist.tsv")
    wordlist_cognate_path = os.path.join(wordlist_cognate_dir,  name + "_wordlist_cognate.tsv")
    glottolog_tree_path = os.path.join(glottolog_tree_dir, name + "_glottolog.tree")
    bin_msa_path = os.path.join(msa_dir, name + ".bin.phy")
    raxml_dir = os.path.join(results_dir, "raxml", name)
    if not os.path.isdir(raxml_dir):
        os.makedirs(raxml_dir)
    raxml_prefix = os.path.join(raxml_dir, "inference")
    pythia_prefix = os.path.join(pythia_dir, name)
    best_tree_path = raxml_prefix + ".raxml.bestTree"

    pipeline.detect_cognates(wordlist_path, wordlist_cognate_path, redo)
    cd = CategoricalData.from_edictor_tsv(wordlist_cognate_path)
    cd.write_msa(bin_msa_path, "bin")
    tree = cd.get_glottolog_tree()
    tree.write(format = 1, outfile = glottolog_tree_path)
    pipeline.run_raxmlng(bin_msa_path, "BIN+G", raxml_prefix, redo)
    pipeline.run_pythia(bin_msa_path, pythia_prefix, redo)
    print(name)
    print("GQ distance")
    print(util.gq_distance(glottolog_tree_path, best_tree_path))
    print("Pythia difficulty score")
    print(util.get_difficulty(pythia_prefix))


