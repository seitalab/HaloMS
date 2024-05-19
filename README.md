# Mapping Adipocyte Interactome Networks with HALOTAG Enrichment Mass Spectrometry

This repository provides the methodology used in our study to visualize latent variables and protein interaction scores for protein structures using AlphaFold2-based Localcolabfold.

## Preparation

1. **Install Localcolabfold**  
   Install [Localcolabfold](https://github.com/YoshitakaMo/localcolabfold) in your local environment.

2. **Load Target Protein Sequences**  
   Retrieve the target protein sequences from [UniProt](https://www.uniprot.org/).

## Setup

```bash
$ poetry install
```

## Calculation (Localcolabfold)

1. **Predict Protein Interaction Structures**  
   Calculate the structure predictions for protein interactions.
   (`localcolabfold/colab-run-multi-seed0.sh`, `localcolabfold/get_score_multi.py`)

2. **Predict Individual Protein Structures and Extract Latent Variables**  
   Calculate the structure predictions for individual proteins and extract latent variables (256 dimensions) from the model.
   (`localcolabfold/colab-run-mono-seed0.sh`, `localcolabfold/af2emb.py`)

3. **Merge Data**  
   Merge the latent variables of proteins, the structure prediction scores from protein interactions, and the pathway data from KEGG.
   (`localcolabfold/merge_score_af2emb.py`)

\*Note: The calculations in step 1, 2 and 3 are extensive, and we provide the calculated results. (`data`)

## Visualization

1. **Dimensionality Reduction and Visualization**  
   Using UMAP, reduce the dimensionality of the latent variables to two dimensions for interacting proteins with each bait protein used in the Halo-MS method. Visualize these reduced variables together with the interaction scores.
   (Run : `visualization/makefig_adsc.ipynb` / Results : `figures`)

   \*Note: The figure generated here is used in the paper.
