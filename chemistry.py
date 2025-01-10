from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import Draw
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

from typing import List


EXAMPLE_COMPOUNDS = [
    # smiles, substructure fingerprint
    "COC1=CC(C(C2=C3C=NC=C2)=CN(C3=O)C)=CC(OC)=C1CN4CCN(CC4)CCOCCOCCOC5=CC(C6=C(N=CS6)C)=CC=C5CNC([C@@H]7C[C@H](CN7C([C@@H](NC(C8(CC8)F)=O)C(C)(C)C)=O)O)=O",
    "ClC(C=C1)=CC=C1C2(CCCCC2)C3=NC4=CC(N5CCN(C(CCCC(NCCN(CCO6)CC6(CC7)CCN7C(CN(N=C8)C=C8C9=CC%10=C(N9)N=CN=C%10C%11=C(C)C(NC(N%12CC(CC(C)(C)C)(O)C%12)=O)=CC(F)=C%11)=O)=O)=O)CC5)=CC=C4C(NCCN)=N3",
    "C(C1=C2C(=CC=C1)N([C@H]1C(=O)NC(=O)CC1)C(=O)N2C)#CCOC1CCN(C[C@H]2CC[C@H](N3C=C(NC(C4=C5N=C(N6C[C@@H]7OC[C@H]6C7)C=CN5N=C4)=O)C(C(F)F)=N3)CC2)CC1",
    "ClC1=CC(O[C@H]2C(C)(C)[C@H](NC(C3=CC=C(C#CC4CCN(C(C[C@H](NC([C@@H]5C[C@@H](O)CN5C([C@@H](C6=CC(C)=NO6)C(C)C)=O)=O)C7=CC=C(C8=C(C)N=CS8)C=C7)=O)CC4)C=C3)=O)C2(C)C)=CC=C1C#N",
]


# input smiles
def get_rdkit_object(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    mol2 = Chem.AddHs(mol)
    return mol2


# NRot bonds
def get_rot_bonds(mol):
    nRot = Lipinski.NumRotatableBonds(mol)
    return nRot


# HBD
def get_HBD(mol):
    hbd = Lipinski.NumHDonors(mol)
    return hbd


# HBA
def get_HBA(mol):
    hba = Lipinski.NumHAcceptors(mol)
    return hba


# Molecular Weight
def get_molecular_weight(mol):
    mw = Chem.Descriptors.MolWt(mol)
    return mw


# TPSA
def get_tpsa(mol):
    tpsa = Chem.rdMolDescriptors.CalcTPSA(mol)
    return tpsa


# conformational search
def conformer_search(mol):
    # Number of conformers to be generated
    numconfs = 25
    max_iter = 50
    mol2 = Chem.addHs(mol)
    params = AllChem.ETKDGv3()
    cids = AllChem.EmbedMultipleConfs(mol2, numconfs, params)
    res = AllChem.MMFFOptimizeMoleculeConfs(mol2)
    return mol2


# wrap everything
def compute_all(smiles):
    obj = get_rdkit_object(smiles)
    nrot = get_rot_bonds(obj)
    hbd = get_HBD(obj)
    hba = get_HBA(obj)
    mw = get_molecular_weight(obj)
    tpsa = get_tpsa(obj)
    props = {"nrot": nrot, "hbd": hbd, "hba": hba, "mw": mw, "tpsa": tpsa}
    return props


compute_all(EXAMPLE_COMPOUNDS[3])
