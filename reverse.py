import os
import pandas as pd

import epitran
from epitran.backoff import Backoff

from ipatok import tokenise

from lingpy import *
from lingpy.sequence.sound_classes import token2class

from tabulate import tabulate
import subprocess

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
                return (epi, "epi")
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
                return (epi, "backoff")
            except Exception as e:
                print("backoff failed")
                print(e)
                return None
    else:
        print(epitran_codes)
        return None


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


def check_ipatok_lexibank_analyzed():
    df = pd.read_csv(os.path.join("resources", "lexibank-analyzed_wordlists", "lexibank-anaylzed_wordlist.tsv"), sep = "\t")
    df.astype("str")
    ipatok_wrong = 0
    lingpy_wrong = 0
    overall = 0
    for i, row in df.iterrows():
        ipa = row["IPA"]
        tokens = row["TOKENS"]
        if ipa == ipa and ipa != "" and ipa != "nan":
            if tokens == tokens and tokens != "" and tokens != "nan":
                overall += 1
                tokens2 = " ".join(tokenise(ipa))
                tokens3 = "".join(ipa2tokens(ipa.replace(" ", "")))
                if tokens2 != tokens:
                    ipatok_wrong += 1
                if tokens3 != tokens:
                    lingpy_wrong += 1
    print("lexibank-analyzed")
    print("error rate ipatok:", str(100 * (ipatok_wrong / overall)))
    print("error rate lingpy:", str(100 * (lingpy_wrong / overall)))


def check_ipatok(df):
    overall = 0
    ipatok_wrong = 0
    lingpy_wrong = 0
    for i, row in df.iterrows():
        ipa = row["rawIPA"]
        tokens = row["IPA"]
        if ipa == ipa and ipa != "" and ipa != "nan":
            if tokens == tokens and tokens != "" and tokens != "nan":
                overall += 1
                tokens = tokens.split(" ")
                tokens2 = tokenise(ipa)
                tokens3 = "".join(ipa2tokens(ipa.replace(" ", "")))
                if tokens2 != tokens:
                    ipatok_wrong += 1
                if tokens3 != tokens:
                    lingpy_wrong += 1
    print("northeuralex")
    print("error rate ipatok:", str(100 * (ipatok_wrong / overall)))
    print("error rate lingpy:", str(100 * (lingpy_wrong / overall)))

def check_epitran(df, relevant_glottocodes):
    results = []
    dolgo_model = rc('dolgo')
    for glottocode in relevant_glottocodes:
        print(glottocode)
        overall = 0
        wrong = 0
        dolgo_wrong = 0
        res = get_epitran(glotto_iso_map[glottocode], epitran_dict)
        assert(res is not None)
        epi = res[0]
        mode = res[1]
        sub_df = df[df["Glottocode"] == glottocode]
        for i, row in sub_df.iterrows():
            ipa = row["rawIPA"]
            form = row["Word_Form"]
            if ipa == ipa and ipa != "" and ipa != "nan" \
                and form == form and form != "" and form != "nan" :
                try:
                    ipa2 = epi.transliterate(form)
                except Exception as e:
                    print(e)
                    wrong += 1
                    continue

                tokens = tokenise(ipa)
                dolgo = [token2class(token, dolgo_model) for token in tokens]
                ipa2 = epi.transliterate(form)
                tokens2 = tokenise(ipa2)
                dolgo2 = [token2class(token, dolgo_model) for token in tokens2]
                overall += 1
                if ipa2 != ipa:
                    wrong +=1
                if dolgo != dolgo2:
                    dolgo_wrong +=1

        results.append([glottocode, "$" + str(round(100 * (wrong / overall), 2)) + " \%$", "$" + str(round(100 * (dolgo_wrong / overall), 2)) + " \%$"])
    print(tabulate(results, tablefmt = "latex_raw", headers=["glottocode", "error rate", "dolgo error rate"]))



df = pd.read_csv(os.path.join("resources", "northeuralex-0.9-forms.tsv"), sep = "\t")
df.astype("str")
glottocodes = list(set(df["Glottocode"]))
epitran_dict = get_epitran_dict()
glotto_iso_map = get_glotto_iso_map()
relevant_glottocodes = []
for glottocode in glottocodes:
    if glottocode not in glotto_iso_map:
        continue
    iso = glotto_iso_map[glottocode]
    if iso in epitran_dict:
        relevant_glottocodes.append(glottocode)

check_ipatok(df)
check_ipatok_lexibank_analyzed()
check_epitran(df, relevant_glottocodes)
