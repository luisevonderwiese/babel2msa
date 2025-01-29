import os
from tabulate import tabulate
from lingpy.compare.sanity import average_coverage
from lingpy import *


def print_coverage(d, ending = ".tsv"):
    print(d)
    res = []
    for f in os.listdir(d):
        if not f.endswith(ending):
            continue
        path = os.path.join(d, f)
        wl = Wordlist(path)
        try:
            res.append([f, average_coverage(wl)])
        except:
            res.append([f, float("nan")])
    print(tabulate(res, tablefmt="pipe", floatfmt=".3f", headers = ["wordlist", "average mutual coverage"]))



print_coverage("results/lexibank-analyzed_families/wordlist_cognate", "_all_wordlist_cognate.tsv")
print_coverage("results/lexibank-analyzed_languagelists/wordlist_cognate", "_all_wordlist_cognate.tsv")
print_coverage("resources/lexibench_wordlists", "_all_wordlist_cognate.tsv")

print_coverage("results/all/wordlist_cognate")
print_coverage("results/main/wordlist_cognate")
print_coverage("results/iecor/wordlist_cognate")


