# Babel2MSA
Project for generating wordlists and character matrices for phylogenetic inference from the Index of [BabelNet](https://babelnet.org/)

## Requierments:
- Install JDK and set `JAVA_HOME`
- Set up the Conda Environment
```
conda env create -f environment.yml
```
- Install `lex_lookup` (for `epitran`) as explained [here](https://github.com/dmort27/epitran)
- Install [BabelNet-API Version 5.3](https://babelnet.org/downloads) in `BabelNet-API-5.3/`
- Place [BabelNet-Index Version 5.0](https://babelnet.org/downloads) in `BabelNet-5.0/`

## Execution:
```
python convert_core_wordnet.py
python experiment.py
python analyze_synsetfilter.py
python completeness_analysis.py
python experiment_lexibench.py
python experiment_lexibank_analyzed.py
python summarize.py
python entropies_lexibench.py
python entropies_lexibank_analyzed.py
python reverse.py
python northeuralex_signal.py
```
