import os
import pandas as pd
import matplotlib.pyplot as plt

ranking_dir = "results/all/ranking"

outdir = "results/analyze_synsetfilter"
if not os.path.isdir(outdir):
    os.makedirs(outdir)

for name in ["filter",  "filter_epitran"]:
    path = os.path.join(ranking_dir, name + ".tsv")
    df = pd.read_csv(path, sep = "\t")
    print(df["count"])
    plt.hist(df["count"], bins = 20)
    plt.xlabel("Number of languages")
    plt.ylabel("Number of synsets")
    plt.yscale("log")
    plt.savefig(os.path.join(outdir, name + "_hist.png"))
    plt.clf()
    plt.close()
