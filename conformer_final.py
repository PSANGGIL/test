from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import Draw
import numpy as np

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

def smi2xyz_split(f_idx, d_smi, split_distance, d_atom, n_conf):
    try:
        split_smi = [s for s in d_smi.split('.')]
        for idx, smiles in enumerate(split_smi):
            smiles.replace('H' , d_atom)
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)
            num_conformers = n_conf
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

        with open(str(f_idx) + '_split_.xyz', 'w') as xyz_file:
            xyz_file.write(f'{len(atom_symbols)}\n')
            xyz_file.write('Generated from SMILES to XYZ conversion (Lowest Energy Conformer)\n')
            for symbol, coord in zip(atom_symbols, lowest_energy_coordinates):
                xyz_file.write(f'{symbol} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
        print(str(f_idx) + 'smi2xyz_split done')
    except:
        print(str(f_idx) + 'smi2xyz_split error')

def smi2xyz(f_idx, smiles, n_conf):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        num_conformers = 100  # Specify the number of conformers to generate
        conformers = AllChem.EmbedMultipleConfs(mol, num_conformers)
        energies = []
        for i, conformer in enumerate(conformers):
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=i)  # Use the index as confId
            ff.Minimize()
            energy = ff.CalcEnergy()
            energies.append(energy)
        lowest_energy_index = energies.index(min(energies))
        conf = mol.GetConformer(lowest_energy_index)

        atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        lowest_energy_coordinates = (conf.GetPositions())

        with open(str(idx) + '_.xyz', 'w') as xyz_file:
            xyz_file.write(f'{len(atom_symbols)}\n')
            xyz_file.write('Generated from SMILES to XYZ conversion (Lowest Energy Conformer)\n')
            for symbol, coord in zip(atom_symbols, lowest_energy_coordinates):
                xyz_file.write(f'{symbol} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
        print(str(f_idx) + 'smi2xyz done')

    except:
        print(str(f_idx) + 'smi2xyz error')
def smi2mol(f_idx, smiles, n_conf):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        num_conformers = n_conf
        conformers = AllChem.EmbedMultipleConfs(mol, num_conformers)
        energies = []
        for i, conformer in enumerate(conformers):
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=i)  # Use the index as confId
            ff.Minimize()
            energy = ff.CalcEnergy()
            energies.append(energy)
        lowest_energy_index = energies.index(min(energies))
        conf = mol.GetConformer(lowest_energy_index)

        #print(energies)
        #print(max(energies))
        #print(min(energies))

        atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        lowest_energy_coordinates = (conf.GetPositions())
        Chem.MolToMolFile(  mol, str(idx) + '_.mol')
        print(str(f_idx) + 'smi2mol done')
    except:
        print(str(f_idx) + 'smi2mol error')

def smi2png(idx, smiles):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    img = Draw.MolToImage(mol)
    img.save( str(idx) + "_.png")
    print('smi2png done')
if __name__ == "__main__":
    # trimol2dimol
    smiles = '[C:14]([C:16](=[O:17])[H:19])([H:20])([H:21])[H:22].[C:1](\[O:4][H:5])([H:7])=[C:13](/[C:15]([H:23])([H:24])[H:25])[H:18].[C:2]([C:3]([H:10])([H:11])[H:12])([H:6])([H:8])[H:9],[C:1]([O:4][H:5])([H:6])([H:7])[C:13]([C:15]([H:23])([H:24])[H:25])([H:18])[H:19].[C:2]([C:3]([H:10])([H:11])[H:12])([H:8])([H:9])[C:16]([C:14]([H:20])([H:21])[H:22])=[O:17]'
    # reactant : tetragonal ring
    smiles = '[C:1](/[C:6](=[C:2]([N:3]([H:8])[H:9])\[N:7]=[C:5](/[N:4]([H:10])[H:11])[H:12])[H:13])([H:14])([H:15])[H:16],[C:1]([C@:6]1([H:13])[C@@:2]2([N:3]([H:8])[H:9])[C@@:5]([N:4]([H:10])[H:11])([H:12])[N:7]21)([H:14])([H:15])[H:16]'
    # product : trigonal ring
    smiles = '[c:1]1([O:5][H:7])[c:3]([H:9])[o:4][c:2]([H:8])[n:6]1,[C@:1]12([O:5][H:7])[C@:2]3([H:8])[O:4][C@@:3]1([H:9])[N:6]23'

    # local(2022.09) done, sever(2023.09) error
    smiles = '[C:1](/[C:5](=[C:6](/[C:2]([H:16])([H:17])[H:18])[C:9](=[C:7]([C:3]([H:19])([H:20])[H:24])[C:4]([H:21])([H:22])[H:23])[H:12])[N:8]([H:10])[H:11])([H:13])([H:14])[H:15],[C:1]([C@:5]1([N:8]([H:10])[H:11])[C@:6]2([C:2]([H:16])([H:17])[H:18])[C:7]([C:3]([H:19])([H:20])[H:24])([C:4]([H:21])([H:22])[H:23])[C@:9]12[H:12])([H:13])([H:14])[H:15]'

    smiles = [s for s in smiles.split(',')]
    n_conf      = 10
    d_atom      = 'F'
    d_distance  = 2

    for idx ,smi in enumerate(smiles):
        smi2png(idx, smi)
        smi2mol(idx, smi, n_conf)
        smi2xyz(idx, smi, n_conf)
        smi2xyz_split(idx, smi, d_distance, d_atom, n_conf)
