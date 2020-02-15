import time, os, re, pyrosetta
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdDepictor, rdMolTransforms
import rdkit.Chem.Draw as Draw
from typing import List, Dict
import os

pyrosetta.init(f'-mute all')

##########################
# TAKEN FROM https://blog.matteoferla.com/2020/02/guess-bond-order-in-rdkit-by-number-of.html
def fix_bond_order(mol: Chem.Mol) -> Chem.Mol:
    """On a Mol where hydrogens are present it guesses bond order."""

    def is_sp2(atom: Chem.Atom) -> bool:
        N_neigh = len(atom.GetBonds())
        symbol = atom.GetSymbol()
        if symbol == 'H':
            return False
        elif symbol == 'N' and N_neigh < 3:
            return True
        elif symbol == 'C' and N_neigh < 4:
            return True
        elif symbol == 'O' and N_neigh < 2:
            return True
        else:
            return False

    def get_other(bond: Chem.Bond, atom: Chem.Atom) -> Chem.Atom:
        """Given an bond and an atom return the other."""
        if bond.GetEndAtomIdx() == atom.GetIdx():  # atom == itself gives false.
            return bond.GetBeginAtom()
        else:
            return bond.GetEndAtom()

    def find_sp2_bonders(atom: Chem.Atom) -> List[Chem.Atom]:
        return [neigh for neigh in find_bonders(atom) if is_sp2(neigh)]

    def find_bonders(atom: Chem.Atom) -> List[Chem.Atom]:
        return [get_other(bond, atom) for bond in atom.GetBonds()]

    def descr(atom: Chem.Atom) -> str:
        return f'{atom.GetSymbol()}{atom.GetIdx()}'

    ## main body of function
    for atom in mol.GetAtoms():
        # print(atom.GetSymbol(), is_sp2(atom), find_sp2_bonders(atom))
        if is_sp2(atom):
            doubles = find_sp2_bonders(atom)
            if len(doubles) == 1:
                # tobedoubled.append([atom.GetIdx(), doubles[0].GetIdx()])
                b = mol.GetBondBetweenAtoms(atom.GetIdx(), doubles[0].GetIdx())
                if b:
                    b.SetBondType(Chem.rdchem.BondType.DOUBLE)
                else:
                    raise ValueError('Issue with:', descr(atom), descr(doubles[0]))
            elif len(doubles) > 1:
                for d in doubles:
                    b = mol.GetBondBetweenAtoms(atom.GetIdx(), d.GetIdx())
                if b:
                    b.SetBondType(Chem.rdchem.BondType.AROMATIC)
                    b.SetIsAromatic(True)
                else:
                    raise ValueError('Issue with:', descr(atom), descr(d))
            elif len(doubles) == 0:
                print(descr(atom), ' is underbonded!')
        else:
            pass
    return mol


#####################################

class Params2SVG:
    glycine = Chem.MolFromSmiles('C(C(=O)O)N')
    AllChem.Compute2DCoords(glycine)
    dimethylglycine = Chem.MolFromSmiles('CC(C)(C(=O))N')
    alanine = Chem.MolFromSmiles('CC(C(=O))N')

    def __init__(self, fullfile):
        self.fullfile = fullfile
        self.name_it()
        if self.name in ('HIP', 'HP2'):
            raise ValueError('causes seg fault')
        print('*' * 10)
        print('*' * 10)
        print(self.name)
        self.make_pdb()
        self.make_mol()
        #print(self.smiles)
        self.make_png()
        self.make_svg()

    def name_it(self):
        n = [line for line in open(self.fullfile) if line.find('NAME') == 0][0]
        self.name = re.match('NAME (\w+)', n).group(1)
        return self

    def make_pdb(self):
        # pyrosetta.rosetta.basic.options.set_file_option('extra_res_fa',os.path.join(path, file))
        #pyrosetta.init(f'-extra_res_fa {self.fullfile} -mute all')
        pose = pyrosetta.rosetta.core.pose.Pose()
        pyrosetta.rosetta.core.pose.make_pose_from_sequence(pose, 'A', 'fa_standard')
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=1, new_res=self.name).apply(pose)
        self.pdbfile = f'pdbs/{self.name}.pdb'
        self.ciffile = f'cifs/{self.name}.cif'
        pose.dump_pdb(self.pdbfile)
        pose.dump_cif(self.ciffile)
        return self

    def make_mol(self):
        try:
            mol = Chem.MolFromPDBFile(self.pdbfile, removeHs=False)
            if mol is None:
                raise ValueError('test_Unreadible.')
            self.mol = fix_bond_order(mol)
            self.mol = Chem.RemoveHs(self.mol)
            AllChem.Compute2DCoords(self.mol)
            self.Lise()
        except BaseException as err:
            print(f'{err.__class__.__name__}: {err}')
            self.mol = Chem.MolFromPDBFile(self.pdbfile)
            self.Lise()
        # glycine as reference for alignment
        AllChem.GenerateDepictionMatching2DStructure(self.mol, self.glycine)
        self.smiles = Chem.MolToSmiles(self.mol)

    def make_svg(self):
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(400, 200)
        rdDepictor.Compute2DCoords(self.mol)
        drawer.DrawMolecule(self.mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        open(f'svgs/{self.name}.svg', 'w').write(svg)

    def make_png(self):
        Draw.MolToFile(self.mol, f'pngs/{self.name}.png')

    def Lise(self):
        """
        Make L amino acid. Modified to account for aldehyde and possible dimethyl on alpha
        TAKEN FROM from https://blog.matteoferla.com/2019/10/rdkit-for-rosetta-plp-ligand-space-as.html
        """
        # query = Chem.MolFromSmiles('CC(C(=O)O)N') #there is no OXT atom.
        if self.mol is None:
            raise ValueError('Mol is None')
        elif self.mol.HasSubstructMatch(self.dimethylglycine):
            query = self.dimethylglycine
            rep = Chem.MolFromSmiles('CC(C)(C(=O)[O-])N')
        elif self.mol.HasSubstructMatch(self.alanine):
            query = self.alanine
            rep = Chem.MolFromSmiles('C[C@@H](C(=O)[O-])N')
        else:
            print('No idea what this is.')
        self.mol = AllChem.ReplaceSubstructs(self.mol, query, rep, replacementConnectionPoint=0)[0]

def main(path) -> List[Params2SVG]:
    # starting = pyrosetta.pose_from_pdb("template.pdb")
    data = []
    for file in os.listdir(path):
        if '.params' not in file:
            continue
        else:
            try:
                fullfile = os.path.join(path, file)
                data.append(Params2SVG(fullfile))
            except BaseException as err:
                print(err.__class__.__name__, err)
    return data

def make_table(data: List[Params2SVG]):
    for model in data:
        if isinstance(model, Params2SVG):
            print(
                f'<tr><td style=>{model.name}</td><td><img src="svgs/{model.name}.svg"/> </td><td>{os.path.split(model.fullfile)[1]}</td><td><code>{model.smiles}</code></td></tr>')
        else:
            err = model["error"].replace("\n", "<br/>")
            print(f'<tr><td> Error </td><td> X </td><td> {model["file"]} </td><td> {err} </td></tr>')

if __name__ == '__main__':
    path = '/Users/matteoferla/rosetta_bin_mac_2019.35.60890_bundle/main/database/chemical/residue_type_sets/fa_standard/residue_types/l-ncaa'
    data = main(path)
    make_table(data)