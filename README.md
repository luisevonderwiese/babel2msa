# Babel2MSA
Project for generating wordlists and character matrices for phylogenetic inference from the Index of [BabelNet](https://babelnet.org/)

## Requirements:
- Install JDK and set `JAVA_HOME`
- Set up the Conda Environment
```
conda env create -f environment.yml
```
- Install `lex_lookup` (for `epitran`) as explained [here](https://github.com/dmort27/epitran)
- Install [BabelNet-API Version 5.3](https://babelnet.org/downloads) in `BabelNet-API-5.3/`
- Place [BabelNet-Index Version 5.0](https://babelnet.org/downloads) in `BabelNet-5.0/` and set `USE_BABELNET_INDICES = True` in `experiment.py`
alternatively:
- Place precompiled data in `results/` and set `USE_BABELNET_INDICES = False` in `experiment.py`

- Download [Glottolog v5.1](https://github.com/glottolog/glottolog) into `resources/glottolog`:
```
cd  resources/
git clone https://github.com/glottolog/glottolog.git
cd glottolog/
git checkout tags/v5.1
```

## Execution:
```
python convert_core_wordnet.py
python experiment.py
# python analyze_synsetfilter.py
python completeness_analysis.py
python summarize.py
python entropies_lexibench.py
python entropies_lexibank_analyzed.py
python reverse.py
python northeuralex_signal.py
```
## Ground Truth Difficulties
Computation of ground truth difficulty scores in [separate repo](https://github.com/luisevonderwiese/difficulty-prediction-training-data/tree/language_data). See there for details.
