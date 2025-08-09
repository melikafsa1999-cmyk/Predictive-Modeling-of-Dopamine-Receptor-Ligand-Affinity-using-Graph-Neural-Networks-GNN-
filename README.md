# Predictive-Modeling-of-Dopamine-Receptor-Ligand-Affinity-using-Graph-Neural-Networks-GNN-


This repository contains the source code and related resources for the study:

Using Machine Learning in the In-Silico Design of Selective Dopamine Receptor Ligands: Advancements in Targeted Therapies for Neurological and Psychiatric Disorders

Our approach integrates edge-conditioned graph neural networks (GNN) with structure-based modeling to predict binding affinities across dopamine receptor subtypes (D1–D5). The workflow includes:

Data preparation and preprocessing (ligand standardization, descriptor generation)

GNN model training and evaluation (with RMSE, R² metrics)

Model interpretability using GNNExplainer

Virtual screening selectivity assessment via Enrichment Factor (EF)

Molecular docking, pharmacophore modeling, and molecular dynamics (MD) simulations for selected ligands

Repository Structure
plaintext
Copy
Edit
├── data/               # Preprocessed datasets (Ki values, SMILES, receptor IDs)
├── models/             # Trained GNN model files (.pt)
├── scripts/            # Python scripts for data processing, training, evaluation
├── results/            # Generated results (plots, EF calculations, interpretability outputs)
├── README.md           # Project description and usage guide
└── requirements.txt    # Python dependencies
Installation
Clone the repository:

bash
Copy
Edit
git clone https://github.com/melikafsa1999-cmyk/Predictive-Modeling-of-Dopamine-Receptor-Ligand-Affinity-using-Graph-Neural-Networks-GNN.git
cd Dopamine-GNN
Install dependencies:

bash
Copy
Edit
pip install -r requirements.txt


Data Availability
The dataset used in this study (BindingDB subset for dopamine receptor ligands) is provided in Data/ and includes:

SMILES strings and their corresponding activity for all ligands




