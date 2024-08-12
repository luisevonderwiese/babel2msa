import os
import pandas as pd

from ipatok import tokenise
from lingpy import *
from categorical import CategoricalData
import util

import matplotlib.pyplot as plt


from it.uniroma1.lcl.babelnet import BabelNet
from it.uniroma1.lcl.babelnet import BabelNetUtils
from it.uniroma1.lcl.babelnet import BabelSynsetID
from it.uniroma1.lcl.jlt.util import UniversalPOS
from it.uniroma1.lcl.babelnet import BabelNetQuery
from it.uniroma1.lcl.jlt.util import Language


def get_languages(language_set):
    assert(language_set in ["all", "iecor", "main"])
    langs = Language.values()
    codes = [util.get_code(lang) for lang in langs]
    iso_glotto_map = util.get_iso_glotto_map()
    if language_set != "all":
        with open(os.path.join("resources", language_set + "_languages.txt"), "r") as lang_file:
            subset_glottocodes = lang_file.read().split("\n")
    selected_langs = []
    for l, code in enumerate(codes):
        if code in iso_glotto_map:
            if language_set == "all" or iso_glotto_map[code][0] in subset_glottocodes:
                selected_langs.append(langs[l])
    result = set(selected_langs)
    return result

def old_filter_synsets(bn, synsetfilter_path, langs, epitran_langs, redo):
    if os.path.isfile(synsetfilter_path) and not redo:
        print("Synset filter present")
        return
    print("Filtering synsets...")
    if len(epitran_langs) == 0:
        print("Set use_epitran=True for generating synset filter")
        assert(False)
    it = bn.getSynsetIterator()
    with open(synsetfilter_path, "w+") as filter_file:
        filter_file.write("\t".join(["babelid","en_concept","form_count", "ipa_count", "ipa_epi_count", "en_sense", "keysense"]) + "\n")
    c = 0
    while it.hasNext():
        synset = it.next()
        if str(synset.getType()) != "Concept":
            continue
        en_sense_opt = synset.getMainSense(Language.EN)
        if en_sense_opt.isPresent():
            en_sense = True
            if en_sense_opt.get().isKeySense():
                keysense = True
            else:
                keysense = False
        else:
            en_sense = False
            keysense = False
        c += 1
        concept = ""
        if en_sense:
            concept = en_sense_opt.get().getSimpleLemma()
        form_count = 0
        ipa_count = 0
        ipa_epi_count = 0
        for l, lang in enumerate(langs):
            sense_opt = synset.getMainSense(lang)
            if not sense_opt.isPresent():
                continue
            form_count += 1
            sense = sense_opt.get()
            if concept == "":
                concept = sense.getSimpleLemma()
            transcriptions = sense.getPronunciations().getTranscriptions()
            if len(transcriptions) != 0:
                ipa_count += 1
            elif lang in epitran_langs:
                ipa_epi_count += 1
        if ipa_count == 0:
            continue
        print(concept)
        with open(synsetfilter_path, "a") as filter_file:
            filter_file.write("\t".join([str(synset.getID()),str(concept),str(form_count), str(ipa_count), str(ipa_epi_count), str(en_sense), str(keysense)]) + "\n")


def filter_synsets(bn, synsetfilter_path):
    if os.path.isfile(synsetfilter_path):
        print("Synset filter present")
        return
    print("Filtering synsets...")
    it = bn.getSynsetIterator()
    langs = Language.values()
    lang_codes = [util.get_code(lang) for lang in langs]
    with open(synsetfilter_path, "w+") as filter_file: 
        filter_file.write("\t".join(["babelid","concept"] + lang_codes) + "\n")
    c = 0
    while it.hasNext():
        synset = it.next()
        if str(synset.getType()) != "Concept":
            continue
        en_sense_opt = synset.getMainSense(Language.EN)
        if en_sense_opt.isPresent():
            en_sense = True
        else:
            en_sense = False
        c += 1
        concept = ""
        if en_sense:
            concept = en_sense_opt.get().getSimpleLemma()
        form_count = 0
        ipa_count = 0
        ipa_epi_count = 0
        result_vector = []
        for l, lang in enumerate(langs):
            sense_opt = synset.getMainSense(lang)
            if not sense_opt.isPresent():
                result_vector.append("00")
                continue
            form_count += 1
            sense = sense_opt.get()
            if concept == "":
                concept = sense.getSimpleLemma()
            transcriptions = sense.getPronunciations().getTranscriptions()
            if len(transcriptions) != 0:
                ipa_count += 1
                result_vector.append("11")
            else:
                result_vector.append("10")
        if ipa_count == 0:
            continue
        print(concept)
        with open(synsetfilter_path, "a") as filter_file:
            filter_file.write("\t".join([str(synset.getID()),str(concept)] + result_vector) + "\n")




def babelids_from_conceptlist(bn, conceptlist_path, babelids_path, langs, epitran_langs, redo):
    if os.path.isfile(babelids_path) and not redo:
        print("Synsets present")
        return
    print("Determining synsets...")
    conceptlist_df = pd.read_csv(conceptlist_path, sep = "\t")
    conceptlist_df = conceptlist_df.astype("str")
    with open(babelids_path, 'w+') as babelids_file:
        babelids_file.write("\t".join(["concept", "babelid", "concepticon_id", "concepticon_gloss"]) + "\n")
    for i, row in conceptlist_df.iterrows():
        concept = row["ENGLISH"].split(" ")[0]
        concept = concept.replace("*", "")
        concept = concept.strip("()")
        print(concept)
        query = BabelNetQuery.Builder(concept).from_(Language.EN).POS(UniversalPOS.valueOf(row["POS"])).to(langs).build()
        synsets = bn.getSynsets(query)
        if len(synsets) == 0 or synsets[0] is None:
            continue
        synset = util.select_synset(synsets, langs, epitran_langs)
        with open(babelids_path, 'a') as babelids_file:
            babelids_file.write("\t".join([concept, str(synset.getID()), row["CONCEPTICON_ID"], row["CONCEPTICON_GLOSS"]]) + "\n")
    
def babelids_from_synsetfilter(bn, synsetfilter_path, babelids_path, num_ids, use_epitran, name, redo):
    if os.path.isfile(babelids_path) and not redo:
        print("Synsets present")
        return
    assert(os.path.isfile(synsetfilter_path))
    print("Retrieving BabelIDs...")
    mode = name.split("_")[1]
    df = pd.read_csv(synsetfilter_path, sep = "\t", dtype={'en_concept':'string'})
    concepticon_map = util.get_concepticon_map()
    if use_epitran:
        df = df.sort_values(['ipa_epi_count'], ascending=[0])
    else:
        df = df.sort_values(['ipa_count'], ascending=[0])
    cnt = 0
    with open(babelids_path, 'w+') as babelids_file:
        babelids_file.write("\t".join(["concept", "babelid", "concepticon_id", "concepticon_gloss"]) + "\n")

    for i, row in df.iterrows():
        if cnt == num_ids:
            break
        if mode == "english" and not row["en_sense"]:
            continue
        if mode == "keysense" and not row["keysense"]:
            continue
        cnt += 1
        concept = row["en_concept"]
        print(concept)
        if concept in concepticon_map:
            concepticon_id, concepticon_gloss = concepticon_map[concept]
        else:
            concepticon_id, concepticon_gloss = ("", "")

        with open(babelids_path, 'a') as babelids_file:
            babelids_file.write("\t".join([concept, row["babelid"], str(concepticon_id), concepticon_gloss]) + "\n")

def ranking_from_synsetfilter(bn, synsetfilter_path, ranking_path, langs, epitran_langs, redo):
    if os.path.isfile(ranking_path) and not redo:
        print("Ranking present")
        return
    assert(os.path.isfile(synsetfilter_path))
    print("Ranking synsets...")
    df = pd.read_csv(synsetfilter_path, sep = "\t")
    df = df.astype("str")
    counts_df =  pd.DataFrame(columns = ["babelid", "concept", "count"])
    lang_codes = [util.get_code(lang) for lang in langs]
    epitran_lang_codes = [util.get_code(lang) for lang in epitran_langs]
    for i, row in df.iterrows():
        count = 0
        for code in lang_codes:
            present = row[code] #00: nothing, 10: form only, 11: form and ipa
            if present.endswith("1"):
                count += 1
            # use epitran and language has epitran and language has form
            elif len(epitran_lang_codes) > 0 and code in epitran_lang_codes and present.startswith("1"):
                count += 1
        counts_df.loc[i] = [row["babelid"], row["concept"], count]
    counts_df = counts_df.sort_values(['count'], ascending=[0])
        
    with open(ranking_path, 'w+') as ranking_file:
        ranking_file.write("\t".join(["babelid", "concept", "count"]) + "\n")

    for i, row in counts_df.iterrows():
        with open(ranking_path, 'a') as ranking_file:
            ranking_file.write("\t".join([row["babelid"], row["concept"], str(row["count"])]) + "\n")




def babelids_from_ranking(bn, ranking_path, babelids_path, num_ids, redo):
    if os.path.isfile(babelids_path) and not redo:
        print("Synsets present")
        return
    assert(os.path.isfile(ranking_path))
    print("Retrieving BabelIDs...")
    concepticon_map = util.get_concepticon_map()
    cnt = 0
    df = pd.read_csv(ranking_path, sep = "\t")
    df = df.astype(str)
    with open(babelids_path, 'w+') as babelids_file:
        babelids_file.write("\t".join(["concept", "babelid", "concepticon_id", "concepticon_gloss"]) + "\n")

    for i, row in df.iterrows():
        if cnt == num_ids:
            break
        cnt += 1
        concept = row["concept"]
        print(concept)
        if concept in concepticon_map:
            concepticon_id, concepticon_gloss = concepticon_map[concept]
        else:
            concepticon_id, concepticon_gloss = ("", "")

        with open(babelids_path, 'a') as babelids_file:
            babelids_file.write("\t".join([concept, row["babelid"], str(concepticon_id), concepticon_gloss]) + "\n")



def languages_statistics(langs, epitran_instances):
    codes = [util.get_code(lang) for lang in langs]
    doculects = util.get_doculects(langs)
    glottocodes = util.get_glottocodes(codes)
    epitran_dict = {}
    with open(os.path.join("resources", "epitran_codes.txt"), "r") as codes_file:
        epitran_codes = codes_file.readlines()
        epitran_codes = [code[:-1] for code in epitran_codes]
    for code in epitran_codes:
        lang = code.split("-")[0]
        if not lang in epitran_dict:
            epitran_dict[lang] = []
        epitran_dict[lang].append(code)

    bitstrings = []
    relevant = []
    for l, lang in enumerate(langs):
        if util.has_code(lang):
            bitstring = "1"
        else:
            bitstring = "0"
        if glottocodes[l] == "":
            bitstring += "0"
        else:
            bitstring += "1"
        if codes[l] in epitran_dict:
            bitstring += "1"
        else:
            bitstring += "0"
        if bitstring == "111":
            relevant.append(glottocodes[l])
        bitstrings.append(bitstring)
    print("Overall", str(len(bitstrings)))
    print("Nothing", str(bitstrings.count("000")))
    print("ISO only", str(bitstrings.count("100")))
    print("Glotto only", str(bitstrings.count("010")))
    #print("Epi only", str(bitstrings.count("001")))
    print("ISO and glotto", str(bitstrings.count("110")))
    print("ISO and epi", str(bitstrings.count("101")))
    #print("Glotto and epi", str(bitstrings.count("011")))
    print("All", str(bitstrings.count("111")))
    print(relevant)



def generate_wordlist(bn, babelids_path, wordlist_path, langs, epitran_instances, redo):
    if os.path.isfile(wordlist_path) and not redo:
        print("Wordlist present")
        return
    print("Generating wordlist...")
    assert(os.path.isfile(babelids_path))
    codes = [util.get_code(lang) for lang in langs]
    doculects = util.get_doculects(langs)
    glottocodes = util.get_glottocodes(codes)
    
    id_df = pd.read_csv(babelids_path, sep = "\t")
    id_df = id_df.astype("str")
    with open(wordlist_path, "w+") as wordlist_file:
        wordlist_file.write("\t".join(["ID","DOCULECT","GLOTTOCODE", "ISO_CODE", "CONCEPT","CONCEPTICON_ID", "CONCEPTICON_GLOSS", "FORM", "IPA", "TOKENS"]) + "\n")

    ID = 0
    for c, row in id_df.iterrows():
        babelid = row["babelid"]
        concept = row["concept"]
        concepticon_id = row["concepticon_id"]
        concepticon_gloss = row["concepticon_gloss"]
        synset = bn.getSynset(BabelSynsetID(babelid))
        form_count = 0
        ipa_count = 0
        ipa_epi_count = 0
        for l, lang in enumerate(langs):
            if glottocodes[l] == "":
                continue
            sense_opt = synset.getMainSense(lang)
            if not sense_opt.isPresent():
                continue
            form_count += 1
            sense = sense_opt.get()
            transcriptions = sense.getPronunciations().getTranscriptions()
            if len(transcriptions) != 0:
                form = str(sense.getSimpleLemma()).replace("_", " ")
                ipa  = str(list(transcriptions)[0].split(",")[0]).strip("[]").strip("/")
                ipa_count += 1
            elif epitran_instances[l] is not None:
                form = str(sense.getSimpleLemma()).replace("_", " ")
                epi = epitran_instances[l]
                try:
                    ipa = epi.transliterate(form)
                    ipa_epi_count += 1
                except:
                    continue
            else:
                continue
            ipa, valid = util.process_ipa(ipa)
            if not valid:
                continue
            tokens = tokenise(ipa, strict=False, replace=True, diphthongs=False, tones=False, unknown=False)
            tokens = list(filter(lambda token: token not in ["ʼ", "`", "´", "", " "], tokens))
            if len(tokens) == 0:
                continue
            tokens = " ".join(tokens)
            print("writing")
            with open(wordlist_path, "a") as wordlist_file:
                wordlist_file.write("\t".join([str(ID), doculects[l], str(glottocodes[l]), codes[l], str(c) + "_" + str(concept), str(concepticon_id), concepticon_gloss, form, ipa, tokens]) + "\n")
            ID += 1
        print(concept, str(form_count), str(ipa_count), str(ipa_epi_count))


def detect_cognates(wordlist_path, wordlist_cognate_path, redo):
    if os.path.isfile(wordlist_cognate_path) and not redo:
        print("Cognates present")
        return
    assert(os.path.isfile(wordlist_path))
    print("Detecting cognates...")
    lex = LexStat(wordlist_path, check = False)
    lex.get_scorer()
    lex.cluster(method="lexstat", threshold=0.6, ref="cognates")
    lex.output('tsv', filename=wordlist_cognate_path.split(".")[0], ignore = "all", prettify = False)

def run_raxmlng(msa_path, model, prefix, redo, args = ""):
    if os.path.isfile(prefix + ".raxml.bestTree") and not redo:
        print("Inferred tree present")
        return
    assert(os.path.isfile(msa_path))
    print("Running inference...")
    command = "./bin/raxml-ng"
    command += " --msa " + msa_path
    command += " --model " + model
    command += " --prefix " + prefix
    command += " --threads auto --seed 2 --force model_lh_impr -blopt nr_safe"
    command += " " + args + " --redo"
    os.system(command)

def run_pythia(msa_path, pythia_prefix, redo):
    if os.path.isfile(pythia_prefix) and not redo:
        print("Pythia results present")
        return
    assert(os.path.isfile(msa_path))
    print("Running pythia...")
    command = "pythia -m " + msa_path + " -o " + pythia_prefix + " -r ./bin/raxml-ng -p predictors/latest.pckl --removeDuplicates -v"
    print(command)
    os.system(command)
    d = util.get_difficulty(pythia_prefix)
    if d != d:
        os.remove(pythia_prefix)
        util.write_padded_msa(msa_path, "temp.phy")
        command = "pythia -m temp.phy -o " + pythia_prefix + " -r ./bin/raxml-ng -p predictors/latest.pckl --removeDuplicates -v"
        print(command)
        os.system(command)
        os.remove("temp.phy")


def sparsity_plot(cd, sparsity_plot_path, redo):
    if os.path.isfile(sparsity_plot_path) and not redo:
        print("Sparsity plot present")
        return
    print("Generating sparsity plot...")
    

    data = [[min(1, len(values)) for values in array] for array in cd.matrix]
    data = [[row[i] for row in data] for i in range(len(data[0]))]
    y = max(4, len(data) / 10)
    factor = len(data) / y 
    x = len(data[0]) / (factor * 10)
    if x < y:
        x = len(data[0]) / factor
    print(str(x), str(y))
    fig, ax = plt.subplots(figsize = (x, y))
    ax.imshow(data, cmap='Greys', interpolation='nearest', aspect = 'auto')
    plt.tight_layout()
    plt.savefig(sparsity_plot_path, bbox_inches='tight') 
    plt.close()

def statistics(wordlist_path, plots_dir, stat_dict, redo):
    if "number of languages" in stat_dict and not redo:
        print("Statistics present")
        return stat_dict
    assert(os.path.isfile(wordlist_path))
    print("Calculating statistics...")

    full_name = stat_dict["name"]
    df = pd.read_csv(wordlist_path, sep = "\t")
    languages = list(set(df["DOCULECT"]))
    concepts = list(set(df["CONCEPT"]))
    stat_dict["number of languages"] = len(languages)
    stat_dict["number of concepts"] = len(concepts)

    concept_nums = []
    for l in languages:
        sub_df = df[df["DOCULECT"] == l]
        num_c = len(list(set(sub_df["CONCEPT"])))
        concept_nums.append(num_c)

    stat_dict["concepts per language"] =  sum(concept_nums) / len(concept_nums)
    
    plt.hist(concept_nums, 20)
    plt.xlabel("Number of concepts")
    plt.ylabel("Number of languages")
    if not os.path.isdir(os.path.join(plots_dir, "num_concepts_hist")):
        os.makedirs(os.path.join(plots_dir, "num_concepts_hist"))
    plt.savefig(os.path.join(plots_dir, "num_concepts_hist", full_name + "_num_concepts_hist.png"))
    plt.close()


    language_nums = []
    for c in concepts:
        sub_df = df[df["CONCEPT"] == c]
        num_l = len(list(set(sub_df["DOCULECT"])))
        language_nums.append(num_l)

    stat_dict["languages per concept"] = sum(language_nums) / len(language_nums)

    plt.hist(language_nums, 20)
    plt.xlabel("Number of languages")
    plt.ylabel("Number of concepts")
    if not os.path.isdir(os.path.join(plots_dir, "num_languages_hist")):
        os.makedirs(os.path.join(plots_dir, "num_languages_hist"))
    plt.savefig(os.path.join(plots_dir, "num_languages_hist", full_name + "_num_languages_hist.png"))
    plt.close()

    return stat_dict



def cognate_statistics(wordlist_cognate_path, plots_dir, stat_dict, redo):
    if "cognate classes per concept" in stat_dict and not redo:
        print("Cogante statistics present")
        return stat_dict
    assert(os.path.isfile(wordlist_cognate_path))
    print("Calculating cognate statistics...")

    full_name = stat_dict["name"]
    df = pd.read_csv(wordlist_cognate_path, sep = "\t")

    class_nums = []
    for c in list(set(df["CONCEPT"])):
        sub_df = df[df["CONCEPT"] == c]
        num_classes = len(list(set(sub_df["COGNATES"])))
        class_nums.append(num_classes)

    stat_dict["cognate classes per concept"] = sum(class_nums) / len(class_nums)

    plt.hist(class_nums, 20)
    plt.xlabel("Number of cognate classes")
    plt.ylabel("Number of concepts")
    if not os.path.isdir(os.path.join(plots_dir, "num_classes_hist")):
        os.makedirs(os.path.join(plots_dir, "num_classes_hist"))
    plt.savefig(os.path.join(plots_dir, "num_classes_hist", full_name + "_num_classes_hist.png"))
    plt.close()

    return stat_dict

