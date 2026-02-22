## Project Overview
The overall focus is on using machine learning and Python tools to analyze biological data, estimate chemical bioactivity, and create reproducible workflows — a key aspect of modern bioinformatics projects.

## Features
### Bioactivity Models
ML models predicting biological activity based on chemical features
Notebooks explore data loading, feature generation, model training, and evaluation

### Interactive Components
Python app(s) (e.g., ML_app.py) which can serve as a backbone for interactive prediction tools

### Visualizations
Pre-generated plots and figures illustrating results and data distributions

### Data & Support
Example chemical bioactivity datasets
Supporting files such as pickled models (.pkl) and SMILES data files

## Contents
Bioactivity_Prediciton_Notebooks/   ← Notebooks for data exploration & modeling  
Mannwhitneyu_Test/                  ← Statistical tests and comparisons  
Plots/                              ← Visual summaries (charts, figures)  
bioactivity_corona_data/            ← Sample dataset (CSV + SMILES)  
ML_app.py                           ← App script  
bioactivity_prediction_app.ipynb    ← Main prediction notebook  
requirements.txt                    ← Python dependencies  
logo.png                            ← Project logo
molecules.smi                       ← Example molecules data

## Setup & Installation
1. Clone the repo
```
git clone https://github.com/saamm/Bioinformatics.git
cd Bioinformatics
```

2. Create a Python environment
```
python3 -m venv venv
source venv/bin/activate
```

3. Install dependencies
```
pip install -r requirements.txt
```

4. Launch notebooks or apps

  Jupyter:
  ```
  jupyter notebook
  ```

  Python app (if configured with Streamlit or similar):
  ```
  streamlit run ML_app.py
  ```

## How to Use This Repo
1. For Educators / Learners
Explore notebooks to see workflows in action
Use example datasets to practice analysis

2. For Developers
Reuse model training functions
Integrate the app script into your own workflows

 3. For Researchers
Adapt bioactivity prediction pipelines
Add your own datasets and compare outcomes

## Acknowledgments
Made with passion for bioinformatics and open science


