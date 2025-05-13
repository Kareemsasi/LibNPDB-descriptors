from flask import Flask, render_template, request
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import QED
from rdkit.Chem import rdMolDescriptors

app = Flask(__name__)

def calculate_descriptors(smiles):
    """
    Calculates chemical descriptors and drug-likeness properties using RDKit.
    Handles invalid SMILES input and returns a dictionary of results.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        try:
            # Calculate descriptors
            molecular_weight = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            num_rings = rdMolDescriptors.CalcNumRings(mol)
            num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

            # Calculate drug-likeness properties
            lipinski_violations = 0
            if Lipinski.NumHAcceptors(mol) > 10:
                lipinski_violations += 1
            if Lipinski.NumHDonors(mol) > 5:
                lipinski_violations += 1
            if logp > 5:
                lipinski_violations += 1
            if molecular_weight > 500:
                lipinski_violations += 1
            ghose_filter = Lipinski.NumHAcceptors(mol) <= 10 and Lipinski.NumHDonors(mol) <= 5 and logp <= 5 and molecular_weight <= 500
            veber_rule = num_rotatable_bonds <= 10 and tpsa <= 100 and rdMolDescriptors.CalcNumAtoms(mol) >= 2
            muegge_rule_pass = Lipinski.NumHAcceptors(mol) <= 10 and Lipinski.NumHDonors(mol) <= 5 and -2 <= logp <= 5 and molecular_weight >= 200 and molecular_weight <= 600
            qed = QED.qed(mol)

            results = {
                "smiles": smiles,
                "physicochemical_properties": {
                    "molecular_weight": molecular_weight,
                    "logP": logp,
                    "TPSA": tpsa,
                    "num_rings": num_rings,
                    "num_rotatable_bonds": num_rotatable_bonds
                },
                "pharmacokinetics": {  #  These are very rough estimates, and ideally should come from a model.
                    "HIA_absorption": "High" if logp < 3 else "Low",  # Very simplified
                    "BBB_penetration": "Yes" if logp < 3 and molecular_weight < 400 else "No", # Simplified
                },
                "drug_likeness": {
                    "Lipinski_violations": lipinski_violations,
                    "Ghose_filter": "Yes" if ghose_filter else "No",
                    "Veber_rule": "Pass" if veber_rule else "Fail",
                    "Muegge_rule_pass": "Yes" if muegge_rule_pass else "No",
                    "QED": qed
                }
            }
            return results
        except Exception as e:
            return {"smiles": smiles, "error": f"Error during calculation: {e}"}
    else:
        return {"smiles": smiles, "error": "Invalid SMILES string"}

@app.route("/", methods=["GET", "POST"])
def index():
    results = None
    if request.method == "POST":
        smiles = request.form["smiles"]
        if smiles:
            results = calculate_descriptors(smiles)
    return render_template("index.html", results=results)

if __name__ == '__main__':
    import os
    port = int(os.environ.get("PORT", 10000))  # Render sets PORT environment variable
    app.run(host='0.0.0.0', port=port)
