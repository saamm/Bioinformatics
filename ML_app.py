import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np



# Molecular descriptor calculator
def desc_calc():

    # Read molecule.smi created earlier
    df = pd.read_csv("molecule.smi", sep="\t", header=None)
    df.columns = ["canonical_smiles", "molecule_chembl_id"]

    fingerprint_list = []
    name_list = []

    for smi, name in zip(df["canonical_smiles"], df["molecule_chembl_id"]):
        mol = Chem.MolFromSmiles(str(smi))

        if mol is None:
            fp = [0] * 881
        else:
            fp_vect = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=881)
            fp = list(fp_vect)

        fingerprint_list.append(fp)
        name_list.append(name)

    # Create dataframe with PaDEL-like format
    columns = [f"PubchemFP{i}" for i in range(881)]
    df_fp = pd.DataFrame(fingerprint_list, columns=columns)
    df_fp.insert(0, "Name", name_list)

    # Save exactly as expected by later code
    df_fp.to_csv("descriptors_output.csv", index=False)

    # Optional cleanup (same as your old PaDEL version)
    os.remove("molecule.smi")

# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# Model building
def build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('corona_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df = pd.concat([molecule_name, prediction_output], axis=1)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)

# Logo image
image = Image.open('logo.png')

st.image(image, use_column_width=True)

# Page title
st.markdown("""
# # Bioactivity Prediction App (Acetylcholinesterase)

This app allows you to predict the bioactivity towards inhibting the `Acetylcholinesterase` enzyme. `Acetylcholinesterase` is a drug target for Alzheimer's disease.

**Credits**
- App built in `Python` + `Streamlit` 
- Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) [[Read the Paper]](https://doi.org/10.1002/jcc.21707).
---
""")

# Sidebar
with st.sidebar.header('1. Upload your CSV data'):
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
    st.sidebar.markdown("""
[Example input file](https://raw.githubusercontent.com/saamm/bioinformatics/main/example_acetylcholinesterase.txt)
""")

if st.sidebar.button('Predict'):
    load_data = pd.read_table(uploaded_file, sep=' ', header=None)
    load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)

    st.header('**Original input data**')
    st.write(load_data)

    with st.spinner("Calculating descriptors..."):
        desc_calc()

    # Read in calculated descriptors and display the dataframe
    st.header('**Calculated molecular descriptors**')
    desc = pd.read_csv('descriptors_output.csv')
    st.write(desc)
    st.write(desc.shape)

    # Read descriptor list used in previously built model
    st.header('**Subset of descriptors from previously built models**')
    Xlist = list(pd.read_csv('descriptor_list.csv').columns)
    desc_subset = desc[Xlist]
    st.write(desc_subset)
    st.write(desc_subset.shape)

    # Apply trained model to make prediction on query compounds
    build_model(desc_subset)
else:
    st.info('Upload input data in the sidebar to start!')
