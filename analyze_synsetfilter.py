import os
import pandas as pd
import matplotlib.pyplot as plt

ranking_dir = "results/all/ranking"

outdir = "results/analyze_synsetfilter"
if not os.path.isdir(outdir):
    os.makedirs(outdir)

for name in ["filter",  "filter_epitran"]:
    print(name)
    path = os.path.join(ranking_dir, name + ".tsv")
    df = pd.read_csv(path, sep = "\t")
    print("Percentage of synsets available in 20 languages or less")
    print((len(df[df["count"] <= 20]) / len(df)) * 100)
    print("Percentage of synsets available in more than 60 languages")
    print((len(df[df["count"] > 60]) / len(df)) * 100)
    print("Largest synset")
    print(max(df["count"]))
    plt.hist(df["count"], bins = 20)
    plt.xlabel("Number of languages")
    plt.ylabel("Number of synsets")
    plt.yscale("log")
    plt.savefig(os.path.join(outdir, name + "_hist.png"))
    plt.clf()
    plt.close()
