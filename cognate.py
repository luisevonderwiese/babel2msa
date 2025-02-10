import csv
import math
import pandas as pd
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


illegal_chars = [",", ";", "[", "]"]

class CognateData:

    def __init__(self, language_ids = [], concept_ids = [], matrix = [], glottocodes = [], glottolog_wrapper = None):
        self.language_ids = language_ids
        self.concept_ids = concept_ids
        #store transposed because easier for converting
        self.matrix = matrix #char_id index, language_id index, values
        self.glottocodes = glottocodes
        self.glottolog_wrapper = glottolog_wrapper


    @classmethod
    def from_edictor_tsv(cls, path, glottolog_wrapper = None):
        df = pd.read_table(path, sep="\t", quoting=csv.QUOTE_NONE)
        df = df[df['DOCULECT'].notna()]
        df = df.astype(str)
        language_ids = list(set(df['DOCULECT']))
        concept_ids = list(set(df['CONCEPT']))
        concept_ids.sort()
        language_ids.sort()
        glottocodes = []
        for language_id in language_ids:
            sub_df = df[df["DOCULECT"] == language_id]
            lang_glottocodes = list(set(sub_df['GLOTTOCODE']))
            if len(lang_glottocodes) == 0:
                glottocodes.append("")
            else: #What to do if > 1 ????
                assert(len(lang_glottocodes) == 1)
                if lang_glottocodes[0] != lang_glottocodes[0]:
                    glottocodes.append("")
                glottocodes.append(lang_glottocodes[0])


        matrix = [[[] for language_idx in range(len(language_ids))] for concept_idx in range(len(concept_ids))]
        for index, row in df.iterrows():
            concept_idx = concept_ids.index(row["CONCEPT"])
            language_idx = language_ids.index(row["DOCULECT"])
            value = row["COGNATES"]
            for illegal_char in illegal_chars:
                value = value.replace(illegal_char,"")
            value_set = matrix[concept_idx][language_idx]
            if value not in value_set:
                value_set.append(value)
        for concept_idx in range(len(concept_ids)):
            for language_idx in range(len(language_ids)):
                matrix[concept_idx][language_idx].sort()
        return cls(language_ids, concept_ids, matrix, glottocodes, glottolog_wrapper)




    def __eq__(self, other):
        if self.language_ids != other.language_ids:
            print("language_ids")
            print(self.language_ids)
            print(other.language_ids)
            return False
        if self.concept_ids != other.concept_ids:
            print("concept_ids")
            print(self.concept_ids)
            print(other.concept_ids)
            return False
        if self.matrix != other.matrix:
            print("matrix")
            for i in range(len(self.matrix)):
                for j in range(len(self.matrix[i])):
                    if self.matrix[i][j] != other.matrix[i][j]:
                        print(self.matrix[i][j])
                        print(other.matrix[i][j])
                        return False
            return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)




    def num_languages(self):
        return len(self.language_ids)

    def num_concepts(self):
        return len(self.concept_ids)


    def get_cognate_classes(self, concept_idx):
        value_set = set()
        entries = self.matrix[concept_idx]
        for entry in entries:
            value_set.update(entry)
        cognate_classes = list(value_set)
        cognate_classes.sort()
        return cognate_classes


    def get_glottolog_tree(self):
        if not self.glottolog_wrapper:
            raise Exception("Glottolog Wrapper needs to be specified")
        return self.glottolog_wrapper.get_tree(self.glottocodes, self.language_ids)

    def encode_bin(self, concept_idx):
        codes = []
        possible_values = self.get_cognate_classes(concept_idx)
        for language_idx in range(self.num_languages()):
            language_values = self.matrix[concept_idx][language_idx]
            if len(language_values) == 0: # missing information
                codes.append("-" * len(possible_values))
                continue
            code = ""
            for value in possible_values: #O(possible_values)
                if value in language_values:
                    code += "1"
                else:
                    code += "0"
            codes.append(code)
        return codes


    def write_bin_msa(self, path):
        sequences = ["" for i in range(self.num_languages())]
        for concept_idx in range(self.num_concepts()): #O(num_concepts * loop_complexity)
            codes = self.encode_bin(concept_idx) #O(num_languages * possible_values)
            if codes == []:
                continue
            for (language_idx, code) in enumerate(codes):
                sequences[language_idx] += code
        if sequences[0] == "":
            return False
        records = [SeqRecord(sequences[language_idx],
                             id=str(self.language_ids[language_idx])) for language_idx in range(self.num_languages())]
        msa = MultipleSeqAlignment(records, annotations={}, column_annotations={})
        with open(path,"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(msa)
        return True


    def bin_column_entropy(self, code):
        code = code.replace("-", "")
        entropy = 0
        for char in ["0", "1"]:
            count = code.count(char)
            if count == 0:
                entropy_x = 0
            else:
                prob = count / len(code)
                entropy_x = prob * math.log2(prob)
            entropy += entropy_x
        return -entropy

    def bin_entropy(self):
        entropy = 0
        for concept_idx in range(self.num_concepts()):
            codes = self.encode_bin(concept_idx)
            for code in codes:
                entropy += self.bin_column_entropy(code)
        return entropy
