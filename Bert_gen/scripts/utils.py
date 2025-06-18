
import pandas as pd 
from tqdm import tqdm
import os
from collections import OrderedDict
from rdkit.Chem.rdchem import Mol

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem  # type: ignore
from rdkit.Chem import rdchem
from rdkit.Chem import rdPartialCharges

class Compound3DKit(object):
    @staticmethod
    def get_MMFF_atom_poses(mol, numConfs=None, return_energy=False):
        try:
            new_mol = Chem.AddHs(mol)
            """
            The following functions are used to generate/optimize 3D conformations of the molecule using MMFF force field.
            """
            mul_conform=AllChem.EmbedMultipleConfs(new_mol, numConfs=numConfs) # type: ignore
            res=AllChem.MMFFOptimizeMoleculeConfs(new_mol) # type: ignore
            index = np.argmin([x[1] for x in res])  
            energy = res[index][1]
            conf= new_mol.GetConformer(id=int(index))
        except Exception as e:
            print(f"[WARNING] 3D conformer optimization failed for molecule. Using 2D coordinates. Error: {e}")
            new_mol = Chem.AddHs(mol)
            AllChem.Compute2DCoords(new_mol) # type: ignore
            energy = 0
            conf = new_mol.GetConformer()
            new_mol.SetProp("_3DStatus", "2D_fallback")  # ← Add tag to molecule
        
        atom_poses = Compound3DKit.get_atom_poses(new_mol, conf)
        if return_energy:
            return new_mol, atom_poses, energy
        else:
            return new_mol, atom_poses

    @staticmethod
    def get_atom_poses(mol, conf): 
        """"
        This function returns the optimal x,y,z using the optimal molecule conformation.
        """ 
        atom_poses = [] 
        for i, atom in enumerate(mol.GetAtoms()):
            if atom.GetAtomicNum() == 0:
                print(f"[WARNING] Atom {i} in molecule {mol.GetProp('_Name')} has atomic number 0. Skipping.")
                return [[0.0, 0.0, 0.0]] * len(mol.GetAtoms())
            pos = conf.GetAtomPosition(i)
            atom_poses.append([pos.x,pos.y,pos.z])
        return atom_poses

def mol_to_geognn_graph_data(mol: Mol, atom_poses, dir_type):
    if len(mol.GetAtoms()) == 0:
        print(f"Invalid molecule: {mol}")
        return None
    data = mol_to_graph_data(mol) 
    data['atom_pos'] = np.array(atom_poses, dtype=np.float32)
    data['bind_length'] = Compound3DKit.get_bond_lengths(data['edges'],data['atom_pos'])
    





def mol_to_geognn_graph_data_MMFF3d(smiles):
    """"
    This function generates several molecule conformations and call mol_to_geognn_graph_data to generate the x,y,z for it using the optimal conformation. 
    """
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if mol is None or mol.GetNumAtoms() == 0:
        print(f"Invalid SMILES: {smiles}")
        return None
    if len(mol.GetAtoms()) <= 400:
        mol,atom_poses=Compound3DKit.get_MMFF_atom_poses(mol, numConfs=10)
    else:
        atom_poses = Compound3DKit.get_2d_atom_poses(mol)
    return mol_to_geognn_graph_data(mol, atom_poses, dir_type='HT') 


        




def pretrainprocess(filename):
    df=pd.read_csv(filename,sep='\t',header=None)
    re = []
    for i,smile in enumerate(df.iloc[:,0]):
        print(f"this is {i}th smile")
        #generate adj matrix for each smile
        atom_list, adj_matrix = mol_to_geognn_graph_data_MMFF3d(smile)
        np.save('data/pretrain/adj_matrix/'+str(i)+'.npy', np.array(adj_matrix))
        re.append([atom_list,'data/adj/'+str(i)+'.npy'])
        r=pd.DataFrame(re,columns=['atom_list','adj_matrix'])
        r.to_csv('data/pretrain/atom_list.csv',index=False)




if __name__== '__main__':
    pretrainprocess()
    
    