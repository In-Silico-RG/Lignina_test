#!pip install rdkit matplotlib

import json
import os
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem  

# List of JSON files
json_files = [
    'LigninStructs_3_fn.json',
    'LigninStructs_4_fn.json',
    'LigninStructs_5_fn.json'
]

output_dir_2d = 'lignin_2d_structures'
output_dir_3d = 'lignin_3d_mol'
os.makedirs(output_dir_2d, exist_ok=True)
os.makedirs(output_dir_3d, exist_ok=True)

for json_file in json_files:
    print(f"\nProcessing {json_file}...")
    with open(json_file, 'r') as f:
        data = json.load(f)
    ligninchains = data.get('ligninchains', [])
    for chain in ligninchains:
        lg_id = chain.get('lg_id', 'unknown')
        smiles = chain.get('smilestring', '')
        # 2D Drawing (SVG, works without Cairo)
        mol2d = Chem.MolFromSmiles(smiles)
        if mol2d:
            drawer = Draw.rdMolDraw2D.MolDraw2DSVG(300, 300)
            drawer.DrawMolecule(mol2d)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            svg_path = os.path.join(output_dir_2d, f"{lg_id}.svg")
            with open(svg_path, "w") as fsvg:
                fsvg.write(svg)
            print(f"Saved 2D SVG for {lg_id} to {svg_path}")
        else:
            print(f"Invalid SMILES for {lg_id}: {smiles}")
        # 3D Structure (MOL)
        mol3d = Chem.MolFromSmiles(smiles)
        if mol3d:
            mol3d = Chem.AddHs(mol3d)
            success = AllChem.EmbedMolecule(mol3d, AllChem.ETKDG())
            if success == 0:
                AllChem.UFFOptimizeMolecule(mol3d)
                mol3d = Chem.RemoveHs(mol3d)
                molfile = os.path.join(output_dir_3d, f"{lg_id}.mol")
                Chem.MolToMolFile(mol3d, molfile)
                print(f"Saved 3D mol file for {lg_id} to {molfile}")
            else:
                print(f"3D embedding failed for {lg_id}: {smiles}")

print(f"\nAll valid 2D SVGs are saved to '{output_dir_2d}', and 3D .mol files to '{output_dir_3d}'.")
