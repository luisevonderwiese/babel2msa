import os
import json
import pandas as pd

for language_set in ["all", "iecor", "lexibank-analyzed_families", "lexibank-analyzed_languagelists", "lexibench", "main"]:
    print(language_set)
    stat_dicts = []

    base_dir = os.path.join("results", language_set)

    statistics_dir = os.path.join(base_dir, "statistics")
    for stats_file in os.listdir(statistics_dir):
        with open(os.path.join(statistics_dir, stats_file)) as json_file:
            stat_dicts.append(json.load(json_file))

    df = pd.DataFrame.from_dict(stat_dicts)
    df = df.sort_values(['name'], ascending=[1])
    with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.expand_frame_repr', False):
        print(df.round(3))

