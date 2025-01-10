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
]


# input smiles
def get_rdkit_object(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    mol2 = Chem.AddHs(mol)
    return mol


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


compute_all(EXAMPLE_COMPOUNDS[0])
