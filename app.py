import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw, AllChem, Crippen, Lipinski
import pandas as pd

st.set_page_config(page_title='Dheepthigah | Molecule Explorer', page_icon='🧬', layout='wide')
st.markdown('''<style>
.stApp{background:radial-gradient(circle at 12% 0%,#132b52 0,#071426 42%,#050d1b 100%);color:#eaf2ff}
[data-testid="stSidebar"]{background:#091a31;border-right:1px solid #24466f}
h1,h2,h3{color:#eaf2ff!important} p,label,.stMarkdown{color:#c8d8ef}
div[data-testid="stMetric"]{background:#102541;border:1px solid #28496e;padding:14px;border-radius:14px}
div.stButton>button,.stDownloadButton>button{background:#1769aa;color:white;border:0;border-radius:10px;min-height:2.7rem}
div.stButton>button:hover,.stDownloadButton>button:hover{background:#2588d2;color:white}
.hero{padding:1.5rem 1.7rem;border:1px solid #28496e;border-radius:18px;background:linear-gradient(120deg,#102b4d,#0a1b32);margin-bottom:1rem}
.eyebrow{color:#66b8ff;letter-spacing:.16em;font-size:.75rem;font-weight:700}
</style>''', unsafe_allow_html=True)

EXAMPLES = {
    'Ethanol':'CCO', 'Phenol':'Oc1ccccc1', 'Benzene':'c1ccccc1',
    'Aspirin':'CC(=O)Oc1ccccc1C(=O)O', 'Acetone':'CC(=O)C',
    'Toluene':'Cc1ccccc1', 'Caffeine':'Cn1c(=O)c2c(ncn2C)n(C)c1=O',
    'Paracetamol':'CC(=O)Nc1ccc(O)cc1', 'Cycloheptane':'C1CCCCCC1',
    '3,5-Dimethylhexane':'CC(C)CC(C)C'
}

def parse_smiles(s):
    return Chem.MolFromSmiles(s) if s.strip() else None

def properties(m):
    return {
        'Molecular formula': rdMolDescriptors.CalcMolFormula(m),
        'Molecular weight (g/mol)': Descriptors.MolWt(m),
        'Exact molecular weight (Da)': Descriptors.ExactMolWt(m),
        'LogP': Crippen.MolLogP(m),
        'TPSA (Å²)': rdMolDescriptors.CalcTPSA(m),
        'H-bond donors': Lipinski.NumHDonors(m),
        'H-bond acceptors': Lipinski.NumHAcceptors(m),
        'Rotatable bonds': Lipinski.NumRotatableBonds(m),
        'Ring count': rdMolDescriptors.CalcNumRings(m),
        'Heavy atoms': m.GetNumHeavyAtoms(), 'Total atoms': m.GetNumAtoms(), 'Bonds': m.GetNumBonds()
    }

def show_3d(mol):
    m = Chem.AddHs(mol)
    p = AllChem.ETKDGv3(); p.randomSeed = 42
    if AllChem.EmbedMolecule(m, p) != 0:
        st.warning('Could not generate a 3D conformer for this structure.'); return
    try: AllChem.MMFFOptimizeMolecule(m, maxIters=200)
    except Exception: pass
    block = Chem.MolToMolBlock(m)
    html = f'''<div id="v" style="width:100%;height:450px;background:#081629;border:1px solid #29496c;border-radius:14px"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script><script>
    const el=document.getElementById('v'); const viewer=$3Dmol.createViewer(el,{{backgroundColor:'#081629'}});
    viewer.addModel(`{block}`,'mol'); viewer.setStyle({{}},{{stick:{{radius:0.16}},sphere:{{scale:0.23}}}});
    viewer.zoomTo();viewer.render();viewer.zoom(0.85);
    </script>'''
    components.html(html, height=470, scrolling=False)

st.markdown('''<div class="hero"><div class="eyebrow">CHEMINFORMATICS LEARNING LAB</div>
<h1 style="margin:.35rem 0">🧬 Molecule Explorer</h1>
<p style="font-size:1.05rem">Explore molecular structures, calculate chemical descriptors, and compare compounds.</p>
<span style="color:#91a9c8">Theme by Dheepthigah · Midnight Blue Edition</span></div>''', unsafe_allow_html=True)
with st.sidebar:
    st.markdown('## Navigation')
    page = st.radio('Choose a section',['Molecule Analyzer','Compound Library','Compare Molecules','About'])
    st.caption('Powered by Python + RDKit')

if page == 'Molecule Analyzer':
    st.subheader('Analyze a molecule')
    left,right=st.columns([1,1])
    with left:
        selected=st.selectbox('Choose an example compound',list(EXAMPLES))
        smiles=st.text_input('Or enter your own SMILES notation',value=EXAMPLES[selected])
        st.button('Analyze molecule',type='primary',use_container_width=True)
        st.caption('SMILES = Simplified Molecular Input Line Entry System.')
    mol=parse_smiles(smiles)
    if mol is None:
        st.error('Invalid or empty SMILES. Please check the notation and try again.')
    else:
        props=properties(mol)
        with right:
            st.markdown('#### 2D molecular structure')
            st.image(Draw.MolToImage(mol,size=(520,360)),use_container_width=True)
            st.caption('Canonical SMILES: '+Chem.MolToSmiles(mol))
        st.markdown('#### Molecular properties')
        items=list(props.items())
        for i in range(0,len(items),4):
            cols=st.columns(4)
            for col,(label,value) in zip(cols,items[i:i+4]):
                col.metric(label,f'{value:.3f}' if isinstance(value,float) else str(value))
        st.markdown('#### Interactive 3D molecular view')
        st.caption('Drag to rotate · Scroll or pinch to zoom. The viewer requires an internet connection.')
        show_3d(mol)
        csv=pd.DataFrame([{'Property':k,'Value':v} for k,v in props.items()]).to_csv(index=False).encode()
        st.download_button('Download properties as CSV',csv,'molecule_properties.csv','text/csv')

elif page == 'Compound Library':
    st.subheader('Starter compound library')
    for name,smiles in EXAMPLES.items():
        with st.expander(name):
            mol=parse_smiles(smiles)
            st.code(smiles)
            st.write(properties(mol))

elif page == 'Compare Molecules':
    st.subheader('Compare compounds')
    chosen=st.multiselect('Choose compounds to compare',list(EXAMPLES),default=['Ethanol','Phenol','Benzene','Aspirin'])
    if chosen:
        rows=[]
        for name in chosen:
            row={'Compound':name,'SMILES':EXAMPLES[name]}; row.update(properties(parse_smiles(EXAMPLES[name]))); rows.append(row)
        df=pd.DataFrame(rows).set_index('Compound')
        st.dataframe(df,use_container_width=True)
        st.bar_chart(df[['Molecular weight (g/mol)']])
        st.download_button('Download comparison CSV',df.to_csv().encode(),'molecule_comparison.csv','text/csv')
    else: st.info('Choose at least one compound.')

else:
    st.subheader('About this project')
    st.markdown('''Molecule Explorer is an educational cheminformatics application built with **Streamlit** and **RDKit**.

- Enter and validate SMILES notation.
- Visualize 2D structures and interactive 3D conformers.
- Calculate formula, molecular weight, exact mass, LogP, TPSA, H-bond donors/acceptors, rotatable bonds, ring count, atom and bond counts.
- Compare example compounds and download CSV reports.

Calculated descriptors are computational estimates; they do not independently establish toxicity, safety, biological activity, or experimental properties.''')
    st.write('Example compounds: '+', '.join(EXAMPLES))

st.markdown('---')
st.markdown('<div style="text-align:center;color:#7892b3;font-size:.8rem">MOLECULE EXPLORER · STUDENT PROJECT · DHEEPTHIGAH</div>',unsafe_allow_html=True)
