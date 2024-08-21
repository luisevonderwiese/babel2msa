import os
import pandas as pd
import iso639
import json
import epitran
from epitran.backoff import Backoff
import subprocess



from it.uniroma1.lcl.babelnet import BabelNet
from it.uniroma1.lcl.babelnet import BabelNetUtils
from it.uniroma1.lcl.babelnet import BabelSynsetID
from it.uniroma1.lcl.jlt.util import UniversalPOS
from it.uniroma1.lcl.babelnet import BabelNetQuery
from it.uniroma1.lcl.jlt.util import Language




def my_english_g2p(self, text):
    text = self.normalize(text).lower()
    try:
        arpa_text = subprocess.check_output(['./bin/lex_lookup', text])
        arpa_text = arpa_text.decode('utf-8')
    except OSError:
        logger.warning('lex_lookup (from flite) is not installed.')
        arpa_text = ''
    except subprocess.CalledProcessError:
        logger.warning('Non-zero exit status from lex_lookup.')
        arpa_text = ''
    # Split on newlines and take the first element (in case lex_lookup
    # returns multiple lines).
    arpa_text = arpa_text.splitlines()[0]
    return self.arpa_to_ipa(arpa_text)

#monkeypatch because I cannot install lex_lookup on server
epitran.flite.FliteLexLookup.english_g2p = my_english_g2p

def get_code(l):
    l_string = str(l).lower()
    try:
        iso_lang = iso639.Lang(l_string)
        pt3 = iso_lang.pt3
        if len(pt3) == 0:
            return l_string
        else:
            return pt3
    except:
        return l_string

def has_code(l):
    l_string = str(l).lower()
    try:
        iso_lang = iso639.Lang(l_string)
        pt3 = iso_lang.pt3
        if len(pt3) == 0:
            return False
        else:
            return True
    except:
        return False


def get_concepticon_map():
    with open(os.path.join("resources", "conceptset.json"), "r") as conceptset_file:
        j = json.load(conceptset_file)
    concepticon_map = {}
    for label in  j["conceptset_labels"]:
        ID = int(j["conceptset_labels"][label][0])
        gloss = j["conceptset_labels"][label][1]
        concepticon_map[label] = (ID, gloss)
    for label in  j["alternative_labels"]:
        ID = int(j["alternative_labels"][label][0])
        gloss = j["alternative_labels"][label][1]
        concepticon_map[label] = (ID, gloss)
    return concepticon_map

def get_iso_glotto_map():
    df = pd.read_csv(os.path.join("resources", "languages.csv"))
    iso_glotto_map = {}
    for i, row in df.iterrows():
        iso = row["ISO639P3code"]
        if iso == iso:
            if not iso in iso_glotto_map:
                iso_glotto_map[iso] = []
            iso_glotto_map[iso].append(row["ID"])
    return iso_glotto_map

def get_glotto_iso_map():
    df = pd.read_csv(os.path.join("resources", "languages.csv"))
    glotto_iso_map = {}
    for i, row in df.iterrows():
        iso = row["ISO639P3code"]
        glotto = row["ID"]
        if iso == iso:
            assert(not glotto in glotto_iso_map)
            glotto_iso_map[glotto] =  iso
    return glotto_iso_map


def get_glottocodes(codes):
    all_glottocodes = []
    iso_glotto_map = get_iso_glotto_map()
    for code in codes:
        if code in iso_glotto_map:
            glottocodes = iso_glotto_map[code]
            assert(len(glottocodes) == 1)
            if glottocodes[0] != glottocodes[0]:
                all_glottocodes.append("")
            else:
                all_glottocodes.append(glottocodes[0])
        else:
            all_glottocodes.append("")
    return all_glottocodes

def get_doculects(langs):
    doculects = []
    for lang in langs:
        doculect = lang.getName()
        doculect = doculect.replace("(", "")
        doculect = doculect.replace(")", "")
        doculect = doculect.replace(" ", "_")
        doculects.append(str(doculect))
    return doculects



def get_epitran_instances(langs):
    print("Loading epitran")
    epitran_instances = []
    epitran_dict = {}
    with open(os.path.join("resources", "epitran_codes.txt"), "r") as codes_file:
        epitran_codes = codes_file.readlines()
        epitran_codes = [code[:-1] for code in epitran_codes]
    for code in epitran_codes:
        lang = code.split("-")[0]
        if not lang in epitran_dict:
            epitran_dict[lang] = []
        epitran_dict[lang].append(code)
    for lang in langs:
        code = get_code(lang)
        if code in epitran_dict:
            epitran_codes = epitran_dict[code]
            print(lang.getName())
            if len(epitran_codes) == 1:
                try:
                    epitran_instances.append(epitran.Epitran(epitran_codes[0]))
                except Exception as e:
                    epitran_instances.append(None)
            else:
                try:
                    if codes[0].startswith("cmn"):
                        epitran_instances.append(Backoff(epitran_codes, cedict_file='resources/cedict_1_0_ts_utf-8_mdbg.txt'))
                    else:
                        epitran_instances.append(Backoff(epitran_codes))
                except Exception as e:
                    epitran_instances.append(None)
        else:
            epitran_instances.append(None)
    return epitran_instances


def process_ipa(ipa):
    phonetic_alphabet = ['a', 'ᴀ', 'ã', 'ɑ', 'á', 'à', 'ā', 'ǎ', 'ụ', 'ū', 'â', 'ɛ', 'æ', 'ɜ', \
            'ɐ', 'ʌ', 'e', 'ᴇ', 'ə', 'ẽ', 'ɘ', 'ɤ', 'è', 'é', 'ē', 'ě', 'ê', 'ɚ', 'Œ', 'ɒ',\
            'œ', 'ɞ', 'ɔ', 'ø', 'ɵ', 'o', 'õ', 'ó', 'ò', 'ō', 'ô', 'y', 'ṳ', 'ʏ', 'ʉ', 'u', \
            'ᴜ', 'ʊ', 'i', 'ɪ', 'ɨ', 'ɿ', 'ʅ', 'ɯ', 'ĩ', 'í', 'ǐ', 'ì', 'î', 'ī', 'ɶ', 'ɷ', \
            'ı', 'ǝ', 'ǒ', 'ĭ', 'ŏ', 'ẽ', 'ä', 'ö', 'ǒ', 'ĭ', 'ŏ', 'ŭ', 'ў', 'ă', 'ĕ', 'ü', \
            'ú', 'ũ', 'ṵ', 'ʮ', 'ɩ', 'ỹ', 'ε', 'ù', 'е', 'ï', 'ǔ', 'ạ', 'ụ', 'ọ', 'ỳ', 'ȯ', \
            'û', 'а', 'ę', 'û', 'ị', 'ý', 'å', 'ǫ', 'ë', 'ạ', 'ḭ', 'ḛ', 'ๅ', 'ố', 'ў', 'ȇ', \
            'ȗ', 'ᴇ', 'ε', 'ṍ', 'ṹ', 'ṳ', 'ŷ', 'ʯ', 't͡s', 't͜s', 'd͡z', 'd͜z', 'ʦ', 'ʣ', 't͡ɕ', \
            't͜ɕ', 'd͡ʑ', 'd͜ʑ', 'ʨ', 'ʥ', 't͡ʃ', 'ʄ', 't͜ʃ', 'd͡ʒ', 'd͜ʒ', 'ʧ', 'ʤ', 'c', 'ɟ', \
            't͡ʂ', 't͜ʂ', 'd͡ʐ', 'd͜ʐ', 'č', 't͡θ', 't͜θ', 'k', 'g', 'q', 'ɢ', 'ɡ', 'x', 'ɣ', \
            'χ', 'ǰ', 'ĵ', 'ḳ', 'ǥ', 'ǵ', 'ḡ', 't͡ʂ', 'Ɉ', 'ʈʂ', 'ɖʐ', 'ʈʂʰ', 'tɕ', 'tɕʰ', \
            'dʑ', 'ts', 'dz', 'tsʰ', 'ǃ', 'ǂ', 'ǁ', 'ǀ', 'ʘ', 'gǃ', 'gǂ', 'gǁ', 'gǀ', 'gʘ',\
            'ǃŋ', 'ǂŋ', 'ǁŋ', 'ǀŋ', 'ʘŋ', '!', '|', 'g!', 'g|', '!ŋ', '|ŋ', 'ɠ', 'ʛ', 'ɸ', \
            'β', 'f', 'p͡f', 'p͜f', 'ƀ', 'p', 'b', 'ɓ', 'р', 'ᵐb', 'ᵐp', 'ḇ', 'bv', 'b͡v', 'pf',\
            'ʔ', 'ħ', 'ʕ', 'h', 'ɦ', 'ḥ', 'Ɂ', 'ʡ', 'ʷ', 'j', 'ɥ', 'ɰ', 'm', 'ɱ', 'ʍ', 'ṃ', \
            'n', 'ȵ', 'ɳ', 'ŋ', 'ɴ', 'ň', 'ń', 'ɲ', '∼', 'ṇ', 'ñ', 'ῃ', 'ņ', 'ṋ', 'ɴ', 'ᶇ', \
            's', 'z', '∫', 'ʃ', 'ʒ', 'ʂ', 'ʐ', 'ç', 'ʝ', 'š', 'ž', 'ɕ', 'ʑ', 'ɧ', 'ś', 'ṣ', \
            'ß', 'ŝ', 'ż', 'ẓ', 'ᶊ', 'ʆ', 'ᶎ', 'ɹ', 'ɻ', 'ʀ', 'ɐ̯', 'ɾ', 'r', 'ʁ', 'ɽ', 'l', \
            'ȴ', 'l', 'ɭ', 'ʎ', 'ʟ', 'ɬ', 'ɮ', 'ɫ', 'ł', 'ɺ', 'ḷ', 'ṛ́', 'ṛ', 'ļ', 'ᵲ', 'ř', \
            'ȓ', 'ṙ', 'ᶉ', 'ᶅ', 't', 'd', 'ȶ', 'ȡ', 'ɗ', 'ʈ', 'ɖ', 'θ', 'ð', 'ŧ', 'þ', 'đ', \
            'т', 'ṱ', 'ṭ', 'ḍ', 'ḏ', 'ţ', 'Ɵ', 'ᶁ', 'ƫ', 'w', 'ʋ', 'v', 'ʙ', 'ⱱ', 'ṿ', 'ṽ', \
            'υ', '+', '₁₁', '₂₂', '¹¹', '²²', '₁₂', '₁₃', '₁₄', '₁₅', '₂₃', '₂₄', '₂₅', '₃₄',\
            '₃₅', '₄₅', '¹²', '¹³', '¹⁴', '¹⁵', '²³', '²⁴', '²⁵', '³⁵', '³⁴', '⁴⁵', '₅₁', '₅₂',\
            '₅₃', '₅₄', '₄₃', '₄₂', '₄₁', '₃₂', '₃₁', '₂₁', '⁵¹', '⁵²', '⁵³', '⁵⁴', '⁴¹', '⁴²', \
            '⁴³', '³¹', '³²', '²¹', '₃₃', '³³', '₄₄', '₅₅', '⁵⁵', '⁴⁴', '⁰', '¹', '²', '³', '⁴', '⁵',\
            '⁻', '₁', '₂', '₃', '₄', '₅', '₆', '₀', '˥', '˦', '˨', '˧', '_', '\#', '◦', '·']
    ipa = ipa.replace('"', '')
    ipa = ipa.replace("ә", "ə")
    ipa = ipa.replace("ş", "ʂ")
    ipa = ipa.replace("ỵ", "ɣ")
    ipa = ipa.replace("ʷ", "w")
    ipa = ipa.replace(" ̝", "")
    ipa = ipa.replace("ˤː", "ʕ")
    ipa = ipa.replace("я", "ʁ")
    ipa = ipa.replace("ǀ", "")
    ipa = ipa.replace("і", "i")
    ipa = ipa.replace("ӏ", "ɪ")
    ipa = ipa.replace("ṯ", "ʈ")
    ipa = ipa.replace("ہ", "")
    ipa = ipa.replace("!", "")
    for char in ipa:
        if char not in phonetic_alphabet:
            return ipa, False
    return ipa, True



def select_synset(synsets, langs, epitran_langs):
    selected_synsets = []
    for synset in synsets:
        if str(synset.getType()) != "Concept":
            continue
        if not synset.getMainSense(Language.EN).isPresent():
            continue
        en_sense = synset.getMainSense(Language.EN).get()
        if en_sense.isKeySense():
            selected_synsets.append(synset)
    if len(selected_synsets) > 0:
        if len(selected_synsets) == 1:
            return selected_synsets[0]
        else:
            synsets = selected_synsets
    synset_lang_counts = dict()
    for synset in synsets:
        count = 0
        for lang in langs:
            if synset.getMainSense(lang).isPresent():
                if lang in epitran_langs:
                    count += 1 
                    continue
                transcriptions = synset.getMainSense(lang).get().getPronunciations().getTranscriptions()
                if len(transcriptions) > 0:
                    count += 1
        synset_lang_counts[synset] = count
    synset_lang_counts = {k: v for k, v in sorted(synset_lang_counts.items(), key=lambda item: item[1], reverse = True)}
    for synset, count in synset_lang_counts.items():
        return synset

    
def get_concept_label(synset, langs):
    en_sense_opt = synset.getMainSense(Language.EN)
    if en_sense_opt.isPresent():
        return en_sense_opt.get().getSimpleLemma()
    else:
        for lang in langs:
            sense_opt = synset.getMainSense(lang)
            if sense_opt.isPresent():
                return sense_opt.get().getSimpleLemma()





def get_difficulty(pythia_prefix):
    with open(pythia_prefix, "r") as outfile:
        lines = outfile.readlines()
        if len(lines) == 0:
            return float("nan")
        return float(lines[0])



def gq_distance(tree_name1, tree_name2):
    if tree_name1 is None or tree_name2 is None:
        return float('nan')
    if tree_name1 != tree_name1 or tree_name2 != tree_name2:
        return float("nan")
    os.system("./bin/qdist " + tree_name1 + " " + tree_name2 + " >out.txt")
    lines = open("out.txt").readlines()
    if len(lines) < 2: #error occurred
        print(lines)
        return float('nan')
    res_q = float(lines[1].split("\t")[-3])
    qdist = 1 - res_q
    os.remove("out.txt")
    return qdist



def write_padded_msa(msa_path, outpath):
    with open(msa_path, "r") as msa_file:
        msa_string = msa_file.read()
    parts = msa_string.split("\n\n")
    lines = parts[-1].split("\n")
    block_size = len(lines[1].split(" ")[-1])
    if block_size == 10:
        padding_size = 10
        append_string = " ----------"
    else:
        padding_size = 10 - block_size
        append_string = "-" * padding_size
    if len(parts) != 1:
        msa_string = "\n\n".join(parts[:-1] + ["\n".join([line + append_string for line in lines[:-1]] + [lines[-1]])])
    else:
        msa_string = "\n".join([lines[0]] + [line + append_string for line in lines[1:-1]] + [lines[-1]])

    parts = msa_string.split("\n")
    sub_parts = parts[0].split(" ")

    msa_string = "\n".join([" ".join(sub_parts[:-1] + [str(int(sub_parts[-1]) + padding_size)])] + parts[1:])

    with open(outpath, "w+") as new_msa_file:
        new_msa_file.write(msa_string)
