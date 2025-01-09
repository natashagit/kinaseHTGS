from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import Draw
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.Chem import AllChem

from typing import List


EXAMPLE_COMPOUNDS = [
    # smiles, substructure fingerprint
    "COC1=CC(C(C2=C3C=NC=C2)=CN(C3=O)C)=CC(OC)=C1CN4CCN(CC4)CCOCCOCCOC5=CC(C6=C(N=CS6)C)=CC=C5CNC([C@@H]7C[C@H](CN7C([C@@H](NC(C8(CC8)F)=O)C(C)(C)C)=O)O)=O",
]


# input smiles
def get_rdkit_object(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    print(mol)
    return mol


# NRot bonds
def get_rot_bonds(mol):
    nRot = Lipinski.NumRotatableBonds(mol)
    print(nRot)
    return nRot


# HBD
def get_HBD(mol):
    hbd = Lipinski.NumHDonors(mol)
    print(hbd)
    return hbd


# HBA
def get_HBA(mol):
    hba = Lipinski.NumHAcceptors(mol)
    print(hba)
    return hba


# conformational search
def conformer_search(mol):
    # Number of conformers to be generated
    numconfs = 25
    max_iter = 50
    mol2 = Chem.addHs(mol)
    params = AllChem.ETKDGv3()
    cids = AllChem.EmbedMultipleConfs(mol2, numconfs, params)
    print(len(cids))
    res = AllChem.MMFFOptimizeMoleculeConfs(mol2)
    return mol2


# wrap everything
def compute_all(smiles):
    obj = get_rdkit_object(smiles)
    nrot = get_rot_bonds(obj)
    hbd = get_HBD(obj)
    hba = get_HBA(obj)
    props = {"nrot": nrot, "hbd": hbd, "hba": hba}
    print(props)
    return props


compute_all(EXAMPLE_COMPOUNDS[0])
