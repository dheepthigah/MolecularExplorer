# Dheepthigah — Molecule Explorer

A Streamlit cheminformatics learning app using RDKit.

## Features
- SMILES input and validation
- Ten example molecules: ethanol, phenol, benzene, aspirin, acetone, toluene, caffeine, paracetamol, cycloheptane, and 3,5-dimethylhexane
- 2D molecular structure rendering and interactive 3D conformer viewer
- Molecular formula, molecular weight, exact mass, LogP, TPSA, H-bond donors/acceptors, rotatable bonds, ring count, atom and bond counts
- Compound comparison and CSV downloads

## Run locally
```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## Deploy on Streamlit Community Cloud
1. Upload `app.py`, `requirements.txt`, and `README.md` to the root of a GitHub repository.
2. Visit https://share.streamlit.io/ and sign in with GitHub.
3. Choose **Create app**, select your repository and branch, and set the main file path to `app.py`.
4. Click **Deploy**.

The 3D viewer loads 3Dmol.js from its CDN and therefore needs an internet connection. RDKit calculations and 2D drawings run in Python.
