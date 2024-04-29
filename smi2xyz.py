from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
from tqdm import tqdm
import numpy as np
import pandas as pd
import re
import pickle
def Mol_split(split_distance, d_smi):
    split_smi = [s for s in d_smi.split('.')]
    for idx, smiles in enumerate(split_smi):
        mol = Chem.MolFromSmiles(smiles)
        num_conformers = 10  # Specify the number of conformers to generate
        conformers = AllChem.EmbedMultipleConfs(mol, num_conformers)
        energies = []
        for i, conformer in enumerate(conformers):
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=i)  # Use the index as confId
            ff.Minimize()
            energy = ff.CalcEnergy()
            energies.append(energy)
        lowest_energy_index = energies.index(min(energies))
        conf = mol.GetConformer(lowest_energy_index)
        translation_vector = idx * split_distance * np.array([1, 1, 1])

        transformation_matrix = np.eye(4)
        transformation_matrix[:3, 3] = translation_vector

        rdMolTransforms.TransformConformer(conf, transformation_matrix)

        if idx == 0:
            atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
            lowest_energy_coordinates = (conf.GetPositions())

        else:
            atom_symbols.extend([atom.GetSymbol() for atom in mol.GetAtoms()])
            lowest_energy_coordinates = np.vstack([lowest_energy_coordinates, (conf.GetPositions())])

    return np.array(atom_symbols), lowest_energy_coordinates

def Map_idx(smi, dummy_atom):
    smi_ = H2DUMMY(smi, dummy_atom)
    qmol    = Chem.MolFromSmiles(smi_)
    ind_map = {}
    for atom in qmol.GetAtoms() :
        map_num = atom.GetAtomMapNum()
        if map_num:
            ind_map[map_num-1] = atom.GetIdx()

    return [ind_map[x] for x in sorted(ind_map)]

def H2DUMMY(smi, dummy_atom):
    return  smi.replace('H', dummy_atom)

def DUMMY2H(arr, dummy_atom):
    arr[arr == dummy_atom] = 'H'
    return arr

def Smi2xyz(split_distance, max_i, i, save_path, smi, dummy_atom, agent, save_xyzfile = False):
    save_xyzfile = True
    mol_fragment = smi.split('.')
    fragments = []
    for f in range(len(mol_fragment)):
        fragments.append([int(match.group(1)) for match in re.finditer(r':(\d+)', mol_fragment[f])])
    d_smi = (H2DUMMY(smi, dummy_atom))
    atom_symbols    = Mol_split(split_distance, d_smi)[0][Map_idx(smi, dummy_atom)]
    coordinates     = Mol_split(split_distance, d_smi)[1][Map_idx(smi, dummy_atom)]
    atom_symbols    = DUMMY2H(atom_symbols, dummy_atom)

    if save_xyzfile == True:
        with open(save_path + str(i).zfill(max_i) + '_' + agent + '_' + dummy_atom + '.xyz', 'w') as reactant_xyz_file:
            reactant_xyz_file.write(f'{len(atom_symbols)}\n')
            reactant_xyz_file.write('Dummy atom : ' + dummy_atom + ' / ' + 'Type : ' + agent+'\n')
            for symbol, coord in zip(atom_symbols, coordinates):
                reactant_xyz_file.write(f'{symbol} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
    else:
        pass

    return atom_symbols, coordinates, fragments

def DB(save_path, data, d_atom, rt, split_distance):
    reaction_data  = Raw_Data(split_distance, save_path, data, d_atom, rt)
    with open(rt + '_data.pickle','wb') as fw:
        pickle.dump(reaction_data, fw)


def Raw_Data(split_distance, save_path, data, d_atom, rt):
    error = []
    atom_dict   = {'H':1, 'C':6, 'N':7, 'O':8}
    reaction = data['reaction']
    smiles   = data[rt]
    num_atoms, charges, fragments, positions, rxn, single_fragment = [], [], [], [], [], []
    for i, smiles_ in enumerate(tqdm(smiles)):
        try:
            if smiles_ == smiles.iloc[i]:
                data_ = (Smi2xyz(split_distance, len(smiles), i, save_path, smiles_, d_atom, rt))
                rxn.append(reaction.iloc[i])
                num_atoms.append(len(data_[0]))
                charges.append([atom_dict[i] for i in data_[0]])
                fragments.append(data_[2])
                positions.append(data_[1])
                if len(fragments) ==1:
                    sf = 1
                else:
                    sf = 0
                single_fragment.append(sf)

            else:
                break
#            if len(data_[2]) >= 4:
#                print(smiles.iloc[i])
#                print(reaction.iloc[i])
        except:
            error.append(i)
    #np.save(rt + '.npy', np.array(error))
    #return  {'num_atoms':num_atoms, 'charges':charges, 'fragments':fragments, 'positions':positions, 'rxn':rxn, 'single_fragment':single_fragment}
    return  {'num_atoms':num_atoms, 'charges':charges, 'fragments':fragments, 'positions':positions, 'rxn':rxn}

if __name__ == "__main__":
    save_path = './test/'
    raw_db = pd.read_csv('RGD1CHNO_smiles.csv', index_col = False)
    raw_db = pd.read_csv('db.csv', index_col = False).sort_values(['index'])
    filter_ = np.load('filter.npy')
    raw_db = raw_db.loc[filter_].loc[:1]
    split_distance = 2
    d_atom = 'F'

    (DB(save_path, raw_db, d_atom, 'reactant', split_distance))
    (DB(save_path, raw_db, d_atom, 'product', split_distance))
