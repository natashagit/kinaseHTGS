from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import Draw
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.Chem import AllChem
import py3Dmol

from typing import List


EXAMPLE_COMPOUNDS = [
    # smiles, substructure fingerprint
    "COC1=CC(C(C2=C3C=NC=C2)=CN(C3=O)C)=CC(OC)=C1CN4CCN(CC4)CCOCCOCCOC5=CC(C6=C(N=CS6)C)=CC=C5CNC([C@@H]7C[C@H](CN7C([C@@H](NC(C8(CC8)F)=O)C(C)(C)C)=O)O)=O",
]


# input smiles
def get_rdkit_object(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
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


# get the 3D image
def MolTo3DView(mol, size=(300, 300), style="stick", surface=False, opacity=0.5):
    """Draw molecule in 3D

    Args:
    ----
        mol: rdMol, molecule to show
        size: tuple(int, int), canvas size
        style: str, type of drawing molecule
               style can be 'line', 'stick', 'sphere', 'carton'
        surface, bool, display SAS
        opacity, float, opacity of surface, range 0.0-1.0
    Return:
    ----
        viewer: py3Dmol.view, a class for constructing embedded 3Dmol.js views in ipython notebooks.
    """
    assert style in ("line", "stick", "sphere", "carton")
    mblock = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=size[0], height=size[1])
    viewer.addModel(mblock, "mol")
    viewer.setStyle({style: {}})
    if surface:
        viewer.addSurface(py3Dmol.SAS, {"opacity": opacity})
    viewer.zoomTo()
    return viewer


# wrap everything
def compute_all(smiles):
    obj = get_rdkit_object(smiles)
    nrot = get_rot_bonds(obj)
    hbd = get_HBD(obj)
    hba = get_HBA(obj)
    props = {"nrot": nrot, "hbd": hbd, "hba": hba}
    return props


compute_all(EXAMPLE_COMPOUNDS[0])
