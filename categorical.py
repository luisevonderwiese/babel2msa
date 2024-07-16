import random
import pandas as pd
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from termcolor import colored
import math
import os
import csv

import glottolog


symbols = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G",
         "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X",
         "Y", "Z", "!", "\"", "#", "$", "%", "&", "\'", "(", ")", "*", "+", ",", "/", ":",
         ";", "<", "=", ">", "@", "[", "\\", "]", "^", "_", "{", "|", "}", "~"]

illegal_chars = [",", ";", "[", "]"]

class CategoricalData:

    def __init__(self, taxon_ids = [], char_ids = [], matrix = [], glottocodes = []):
        self.taxon_ids = taxon_ids
        self.char_ids = char_ids
        #store transposed because easier for converting
        self.matrix = matrix #char_id index, taxon_id index, values
        self.glottocodes = glottocodes


    @classmethod
    def from_edictor_tsv(cls, path):
        df = pd.read_csv(path, sep="\t", quoting=csv.QUOTE_NONE)
        df = df[df['DOCULECT'].notna()]
        df = df.astype(str)
        taxon_ids = list(set(df['DOCULECT']))
        char_ids = list(set(df['CONCEPT']))
        char_ids.sort()
        taxon_ids.sort()
        glottocodes = []
        for taxon_id in taxon_ids:
            sub_df = df[df["DOCULECT"] == taxon_id]
            lang_glottocodes = list(set(sub_df['GLOTTOCODE']))
            if len(lang_glottocodes) == 0:
                glottocodes.append("")
            else: #What to do if > 1 ????
                assert(len(lang_glottocodes) == 1)
                if lang_glottocodes[0] != lang_glottocodes[0]:
                    glottocodes.append("")
                glottocodes.append(lang_glottocodes[0])

        matrix = [[[] for taxon_idx in range(len(taxon_ids))] for char_idx in range(len(char_ids))]
        for index, row in df.iterrows():
            char_idx = char_ids.index(row["CONCEPT"])
            taxon_idx = taxon_ids.index(row["DOCULECT"])
            value = row["COGNATES"]
            for illegal_char in illegal_chars:
                value = value.replace(illegal_char,"")
            value_set = matrix[char_idx][taxon_idx]
            if value not in value_set:
                value_set.append(value)
        for char_idx in range(len(char_ids)):
            for taxon_idx in range(len(taxon_ids)):
                matrix[char_idx][taxon_idx].sort()
        return cls(taxon_ids, char_ids, matrix, glottocodes)





    def __eq__(self, other):
        if self.taxon_ids != other.taxon_ids:
            print("taxon_ids")
            print(self.taxon_ids)
            print(other.taxon_ids)
            return False
        if self.char_ids != other.char_ids:
            print("char_ids")
            print(self.char_ids)
            print(other.char_ids)
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

    def num_taxa(self):
        return len(self.taxon_ids)

    def num_chars(self):
        return len(self.char_ids)

    def max_values(self):
        return len(self.max_values_counts()) - 1

    def max_values_counts(self):
        counts = []
        for char_idx in range(self.num_chars()):
            num = len(self.get_possible_values(char_idx))
            while(len(counts) < num + 1):
                counts.append(0)
            counts[num] += 1
        return counts

    def num_sites_bin(self):
        max_values_counts = self.max_values_counts()
        print(max_values_counts)
        return max_values_counts[0] + sum([i * max_values_counts[i] for i in range(1, len(max_values_counts))])

    def sites_per_char(self):
        return self.num_sites_bin() / self.num_chars()


    def size(self):
        return self.num_taxa() * self.num_chars()

    def polymorphic_char_ratio(self):
        cnt = 0
        for c in range(self.num_chars()):
            for t in range(self.num_taxa()):
                if len(self.matrix[c][t]) > 1:
                    cnt += 1
                    break
        return cnt / self.num_chars()


    def get_possible_values(self, char_idx):
        value_set = set()
        entries = self.matrix[char_idx]
        for entry in entries:
            value_set.update(entry)
        possible_values = list(value_set)
        possible_values.sort()
        return possible_values


    def bin_cross_validation_data(self, path):
        random.seed(2)
        ratio = 0.6
        num_samples = 10
        if not os.path.isdir(path):
            os.makedirs(path)
        num_chars_train = math.ceil(self.num_chars() * ratio)
        for t in range(num_samples):
            l = [_ for _ in range(self.num_chars())]
            random.shuffle(l)
            train_indices = l[:num_chars_train]
            train_chars = []
            test_chars = []
            train_matrix = []
            test_matrix = []
            for char_idx in range(self.num_chars()):
                if char_idx in train_indices:
                    train_chars.append(self.char_ids[char_idx])
                    train_matrix.append(self.matrix[char_idx])
                else:
                    test_chars.append(self.char_ids[char_idx])
                    test_matrix.append(self.matrix[char_idx])
            train_cd = CategoricalData(self.taxon_ids, train_chars, train_matrix, self.glottocodes)
            train_cd.write_msa(os.path.join(path, "bin_cv_train_" + str(t) + ".phy"), "bin")
            test_cd = CategoricalData(self.taxon_ids, test_chars, test_matrix, self.glottocodes)
            test_cd.write_msa(os.path.join(path, "bin_cv_test_" + str(t) + ".phy"), "bin")

    def get_glottolog_tree(self):
        return glottolog.get_tree(self.glottocodes, self.taxon_ids)

    def split_families(self):
        languages_dict = {}
        glottocodes_dict = {}
        matrix_dict = {}
        chars_dict = {}
        for i, glottocode in enumerate(self.glottocodes):
            if glottocode == "":
                family = "UNKNOWN"
            else:
                family = glottolog.get_family(glottocode)
            if family not in languages_dict:
                languages_dict[family] = []
                glottocodes_dict[family] = []
                matrix_dict[family] = [[] for char_idx in range(self.num_chars())]
            languages_dict[family].append(self.taxon_ids[i])
            glottocodes_dict[family].append(glottocode)
            for char_idx in range(self.num_chars()):
                matrix_dict[family][char_idx].append(self.matrix[char_idx][i])
        chars_dict = {}
        for family, matrix in matrix_dict.items():
            new_matrix = []
            chars_dict[family] = []
            for char_idx in range(self.num_chars()):
                for language_idx in range(len(matrix[char_idx])):
                    if len(matrix[char_idx][language_idx]) > 0:
                        new_matrix.append(matrix[char_idx])
                        chars_dict[family].append(self.char_ids[char_idx])
                        break
            matrix_dict[family] = new_matrix
        result = {}
        for family, glottocodes in glottocodes_dict.items():
            result[family] = CategoricalData(languages_dict[family], chars_dict[family], matrix_dict[family], glottocodes)

        return result




    def encode_bin(self, char_idx): #O(num_taxa * possible_values)
        codes = []
        possible_values = self.get_possible_values(char_idx)
        for taxon_idx in range(self.num_taxa()): #O(num_taxa * loop_complexity)
            taxon_values = self.matrix[char_idx][taxon_idx]
            if len(taxon_values) == 0: # missing information
                codes.append("-" * len(possible_values))
                continue
            code = ""
            for value in possible_values: #O(possible_values)
                if value in taxon_values:
                    code += "1"
                else:
                    code += "0"
            codes.append(code)
        return codes

    def encode_bin_part(self, char_idx, num_values): #O(num_taxa * possible_values)
        if len(self.get_possible_values(char_idx)) != num_values:
            return []
        codes = []
        possible_values = self.get_possible_values(char_idx)
        for taxon_idx in range(self.num_taxa()): #O(num_taxa * loop_complexity)
            taxon_values = self.matrix[char_idx][taxon_idx]
            if len(taxon_values) == 0: # missing information
                codes.append("-" * len(possible_values))
                continue
            code = ""
            for value in possible_values: #O(possible_values)
                if value in taxon_values:
                    code += "1"
                else:
                    code += "0"
            codes.append(code)
        return codes

    def encode_multi(self, char_idx):
        codes = []
        possible_values = self.get_possible_values(char_idx)
        assert(len(possible_values) <= len(symbols))
        symbol_dict = {value: symbol for value, symbol in zip(possible_values, symbols)}
        for taxon_idx in range(self.num_taxa()):
            taxon_values = self.matrix[char_idx][taxon_idx]
            assert(len(taxon_values) <= 1)
            if len(taxon_values) == 0: # missing information
                codes.append("-")
            else:
                codes.append(symbol_dict[taxon_values[0]])
        return codes


    def encode_ambig(self, char_idx, max_values):
        binary_codes = self.encode_bin(char_idx)
        ambig_codes = []
        for binary_code in binary_codes:
            if binary_code.startswith("-"):
                ambig_codes.append("-")
                continue
            diff = max_values - len(binary_code)
            padd = "0" * diff
            binary_code = padd + binary_code
            idx = int(binary_code, 2)
            assert(idx < pow(2, max_values))
            ambig_codes.append(symbols[idx])
        return ambig_codes

    def encode_prototype(self, char_idx, max_values):
        if len(self.get_possible_values(char_idx)) > max_values:
            return []
        binary_codes = self.encode_bin(char_idx)
        ambig_codes = []
        for binary_code in binary_codes:
            if binary_code.startswith("-"):
                ambig_codes.append("-")
                continue
            diff = max_values - len(binary_code)
            assert(diff >= 0)
            padd = "0" * diff
            binary_code = padd + binary_code
            idx = int(binary_code, 2)
            assert(idx > 0)
            idx = idx - 1
            assert(idx < pow(2, max_values))
            ambig_codes.append(symbols[idx])
        return ambig_codes


    def encode_prototype_part(self, char_idx, num_values):
        if len(self.get_possible_values(char_idx)) != num_values:
            return []
        binary_codes = self.encode_bin(char_idx)
        ambig_codes = []
        for binary_code in binary_codes:
            assert(len(binary_code) == num_values)
            if binary_code.startswith("-"):
                ambig_codes.append("-")
                continue
            idx = int(binary_code, 2)
            assert(idx > 0)
            idx = idx - 1
            assert(idx < pow(2, num_values))
            ambig_codes.append(symbols[idx])
        return ambig_codes

    def get_msa(self, msa_type): #O(num_chars * num_taxa * max(possible_values))
        if msa_type == "ambig":
            max_values = self.max_values()
            if pow(2, max_values) > len(symbols):
                print(colored("ambig MSA cannot be created for dataset with " +  str(max_values) + " max_values", "yellow"))
                return None
            if max_values < 2:
                print(colored("ambig MSA cannot be created for dataset with " +  str(max_values) +  " < 2 max_values", "yellow"))
                return None
        if msa_type == "prototype":
            max_values = min(int(math.floor(math.log(len(symbols), 2))), self.max_values())
        if msa_type == "multi":
            max_values = self.max_values()
            if max_values > len(symbols):
                print(colored("multi MSA cannot be created for dataset with " + str(max_values) + " > 64 max_values", "yellow"))
                return None
            if max_values < 2:
                print(colored("multi MSA cannot be created for dataset with " + str(max_values) + " < 2 max_values", "yellow"))
                return None
            if self.is_mutlistate():
                print(colored("multi MSA cannot be created for multistate dataset", "yellow"))
                return None
        if msa_type.startswith("prototype_part") or msa_type.startswith("bin_part"):
            num_values = int(msa_type.split("_")[-1])

        sequences = ["" for i in range(self.num_taxa())]
        for char_idx in range(self.num_chars()): #O(num_chars * loop_complexity)
            if msa_type == "bin":
                codes = self.encode_bin(char_idx) #O(num_taxa * possible_values)
            elif msa_type == "multi":
                codes = self.encode_multi(char_idx)
            elif msa_type == "ambig":
                codes = self.encode_ambig(char_idx, max_values)
            elif msa_type == "prototype":
                codes = self.encode_prototype(char_idx, max_values)
            elif msa_type.startswith("prototype_part"):
                codes = self.encode_prototype_part(char_idx, num_values)
            elif msa_type.startswith("bin_part"):
                codes = self.encode_bin_part(char_idx, num_values)
            else:
                print("MSA type", msa_type, "invalid")
                return None
            if codes == []:
                continue
            for (taxon_idx, code) in enumerate(codes):
                sequences[taxon_idx] += code
        if sequences[0] == "":
            return None
        #if msa_type.startswith("prototype_part") and len(sequences[0]) < self.num_taxa():
        #    return None
        #if msa_type.startswith("bin_part") and len(sequences[0]) < self.num_taxa() * num_values:
        #    return None
        if (msa_type.startswith("prototype_part") or msa_type.startswith("bin_part")) and len(sequences[0]) == 0:
            return None
        records = [SeqRecord(sequences[taxon_idx],
                             id=str(self.taxon_ids[taxon_idx].replace(" ", "_"))) for taxon_idx in range(self.num_taxa())]
        msa = MultipleSeqAlignment(records, annotations={}, column_annotations={})
        return msa

    def write_msa(self, path, msa_type):
        if msa_type == "catg_bin":
            self.write_catg_msa(path)
            return True
        elif msa_type == "catg_multi":
            return self.write_multi_catg_msa(path)
        msa =  self.get_msa(msa_type)
        if msa is None:
            return False
        with open(path,"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(msa)
        return True


    def write_catg_msa(self, path):
        num_sites = sum([len(self.get_possible_values(char_idx)) for char_idx in range(self.num_chars())])
        probs = [[] for _ in range(num_sites)]
        cols = ["" for _ in range(num_sites)]
        i = 0
        for char_idx in range(self.num_chars()):
            possible_values = self.get_possible_values(char_idx)
            for value in possible_values:
                for taxon_idx in range(self.num_taxa()):
                    taxon_values = self.matrix[char_idx][taxon_idx]
                    if value in taxon_values:
                        # only in this case 1 is more probable than 0
                        if len(taxon_values) == 1:
                            cols[i] += "1"
                        else:
                            cols[i] += "0"
                        one_prob = 1 / len(taxon_values)
                    else:
                        cols[i] += "0"
                        # missing data
                        if len(taxon_values) == 0:
                            probs[i].append("1.0,1.0")
                            continue
                            #alternative:
                            #one_prob = 1 / len(possible_values)
                        else:
                            one_prob = 0.0
                    zero_prob = 1.0 - one_prob
                    one_prob = round(one_prob, 3)
                    zero_prob = round(zero_prob, 3)
                    probs[i].append(str(zero_prob)+","+str(one_prob))
                i += 1

        outfile = open(path, "w+")
        outfile.write(str(self.num_taxa()) + " ")
        outfile.write(str(num_sites) + "\n")
        outfile.write(" ".join(self.taxon_ids))
        outfile.write("\n")
        for i in range(num_sites):
            outfile.write(cols[i] + " ")
            outfile.write(" ".join(probs[i]))
            outfile.write("\n")
        outfile.close()


    def write_multi_catg_msa(self, path):
        probs = []
        cols = []
        max_values = self.max_values()
        if max_values > len(symbols):
            print(colored("catg_multi MSA cannot be created for dataset with " + str(max_values) + " > 64 max_values", "yellow"))
            return False
        if max_values < 2:
            print(colored("catg_multi MSA cannot be created for dataset with " + str(max_values) + " < 2 max_values", "yellow"))
            return False
        for char_idx in range(self.num_chars()):
            possible_values = self.get_possible_values(char_idx)
            symbol_dict = {value: symbol for value, symbol in zip(possible_values, symbols)}
            char_sequence = ["" for _ in range(self.num_taxa())]
            char_probs = ["" for _ in range(self.num_taxa())]
            for taxon_idx in range(self.num_taxa()):
                taxon_values = self.matrix[char_idx][taxon_idx]
                if len(taxon_values) == 0:
                    char_sequence[taxon_idx] = "-"
                    prob_vec = ["1.0" for _ in range(max_values)]
                    char_probs[taxon_idx] = (",".join(prob_vec))
                    continue
                prob_vec = ["0.0" for _ in range(max_values)]
                prob = str(round(1 / len(taxon_values), 3))
                for i, value in enumerate(possible_values):
                    if value in taxon_values:
                        if char_sequence[taxon_idx] == "":
                            char_sequence[taxon_idx] = symbol_dict[value]
                        prob_vec[i] = prob
                char_probs[taxon_idx] = (",".join(prob_vec))
            cols.append("".join(char_sequence))
            probs.append(char_probs)

        outfile = open(path, "w+")
        outfile.write(str(self.num_taxa()) + " ")
        outfile.write(str(self.num_chars()) + "\n")
        outfile.write(" ".join(self.taxon_ids))
        outfile.write("\n")
        for i in range(self.num_chars()):
            outfile.write(cols[i] + " ")
            outfile.write(" ".join(probs[i]))
            outfile.write("\n")
        outfile.close()
        return True

    def get_random_sample(self, seed):
        random.seed(seed)
        sample = CategoricalData()
        sample.char_ids = self.char_ids
        sample.taxon_ids = self.taxon_ids
        sample.glottocodes = self.glottocodes
        sample.matrix = [[[] for taxon_idx in range(self.num_taxa())] for char_idx in range(self.num_chars())]
        for char_idx in range(self.num_chars()):
            for taxon_idx in range(self.num_taxa()):
                if len(self.matrix[char_idx][taxon_idx]) != 0:
                    sample.matrix[char_idx][taxon_idx] = [random.choice(self.matrix[char_idx][taxon_idx])]
        return sample

    def is_single_state(self):
        for char_idx in range(self.num_chars()):
            for taxon_idx in range(self.num_taxa()):
                if len(self.matrix[char_idx][taxon_idx]) > 1:
                    return False
        return True



    def get_value_number_counts(self):
        counts = []
        for char_idx in range(self.num_chars()):
            for taxon_idx in range(self.num_taxa()):
                num = len(self.matrix[char_idx][taxon_idx])
                while(len(counts) < num + 1):
                    counts.append(0)
                counts[num] += 1
        return counts

    def get_value_number_matrix(self):
        counts = [[] for i in range(self.max_values() + 1)]
        for char_idx in range(self.num_chars()):
            num_possible_values = len(self.get_possible_values(char_idx))
            for taxon_idx in range(self.num_taxa()):
                num = len(self.matrix[char_idx][taxon_idx])
                while(len(counts[num_possible_values]) < num + 1):
                    counts[num_possible_values].append(0)
                counts[num_possible_values][num] += 1
        return counts


    def get_multistate_ratio(self):
        return sum(self.get_value_number_counts()[2:]) / self.size()

    def is_mutlistate(self):
        return len(self.get_value_number_counts()) > 2
