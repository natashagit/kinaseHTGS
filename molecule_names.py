import json

data = {
    "Molecule1": "COC1=CC(C(C2=C3C=NC=C2)=CN(C3=O)C)=CC(OC)=C1CN4CCN(CC4)CCOCCOCCOC5=CC(C6=C(N=CS6)C)=CC=C5CNC([C@@H]7C[C@H](CN7C([C@@H](NC(C8(CC8)F)=O)C(C)(C)C)=O)O)=O"
}

with open("molecules.json", "w") as f:
    json.dump(data, f)
