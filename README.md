# Display of the preset Rosetta NCAAs

What exactly are the non-canonical amino acids in the Rosetta database folder?

Names like V02 really do not help... So there is only one way to find out, manually generate a figure for each by recycling a large amount of code.

## Script

The file `generate.py` does the work. To do this, it used both pyrosetta and RDKit. In the case of RDKit the bond order is guessed using a function for a [blog post of mine](https://blog.matteoferla.com/2020/02/guess-bond-order-in-rdkit-by-number-of.html).

In the case of pyrosetta, the mover `pyrosetta.rosetta.protocols.simple_moves.MutateResidue` is used as the argument `new_res` can be a custom parameterised residue. Shamefully, I do not know how one can pass the argument `extra_res_fa` to pyrosetta after it gets initialised as setting it as an option afterward simply failed.

A cool detail is that for once I am not using PyMOL. In fact, the line `pyrosetta.rosetta.core.pose.make_pose_from_sequence(pose, 'A', 'fa_standard')` generates the pose from scratch.

## Gallery

| 3-Letter | File | SMILES | Image |
| --- | --- | --- | --- | --- |
| A34 | 2-aminomethyl-phenylalanine.params | `N[C@@H](Cc1ccccc1C[NH3+])C(=O)[O-]` | ![model.name](pngs/A34.png) |
| A33 | 2-amino-heptanoic_acid.params | `CCCCC[C@H](N)C(=O)[O-]` | ![model.name](pngs/A33.png) |
| Error | YPN.params | X | (RuntimeError)  <br/><br/>File: /Volumes/MacintoshHD3/benchmark/W.fujii.release/rosetta.Fujii.release/_commits_/main/source/src/core/conformation/Residue.cc:1365<br/>[ ERROR ] UtilityExitException<br/>ERROR: Unable to fill in missing atoms.<br/><br/> |
| NVL | NVL.params | `CCC[C@H](N)C(=O)[O-]` | ![model.name](pngs/NVL.png) |
| Error | 4J5.params | X | (AttributeError)  'NoneType' object has no attribute 'GetAtoms' |
| C94 | trifluoro-leucine_ent2.params | `CC(C[C@H](N)C(=O)[O-])C(F)(F)F` | ![model.name](pngs/C94.png) |
| A24 | 2-amino-2-phenylbutyric_acid.params | `CC.N[C@H](C(=O)[O-])c1ccccc1` | ![model.name](pngs/A24.png) |
| B27 | 4-methyl-phenylalanine.params | `Cc1ccc(cc1)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/B27.png) |
| A20 | 2-allyl-glycine.params | `C=CC[C@H](N)C(=O)[O-]` | ![model.name](pngs/A20.png) |
| Error | 5-bromo-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| B57 | alpha-methyl-leucine.params | `CC(C)C.C[C@H](N)C(=O)[O-]` | ![model.name](pngs/B57.png) |
| Error | 7-methyl-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| Error | 4.5-dihydroxy-isoleucine.params | X | (AttributeError)  'NoneType' object has no attribute 'GetAtoms' |
| C26 | homocysteine.params | `N[C@@H](CCS)C(=O)[O-]` | ![model.name](pngs/C26.png) |
| A31 | 2-amino-5-phenyl-pentanoic_acid.params | `N[C@@H](CCCc1ccccc1)C(=O)[O-]` | ![model.name](pngs/A31.png) |
| Error | 0TD.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| Error | 4-amino-tetrahydrothiopyran-4-carboxylic_acid.params | X | (AttributeError)  'NoneType' object has no attribute 'GetAtoms' |
| B97 | beta-chloro-alanine.params | `N[C@@H](CCl)C(=O)[O-]` | ![model.name](pngs/B97.png) |
| B58 | alpha-methyl-phenylalanine.params | `C[C@H](N)C(=O)[O-].Cc1ccccc1` | ![model.name](pngs/B58.png) |
| V02 | V02.params | `N[C@H](C(=O)[O-])c1ccc(O)cc1` | ![model.name](pngs/V02.png) |
| B19 | 4-fluoro-proline.params_rot | `CC(F)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/B19.png) |
| C00 | beta-cyclohexyl-alanine.params | `N[C@@H](CC1CCCCC1)C(=O)[O-]` | ![model.name](pngs/C00.png) |
| Error | 4-fluoro-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| A80 | 3-hydroxy-tyrosine.params | `N[C@@H](Cc1ccc(O)c(O)c1)C(=O)[O-]` | ![model.name](pngs/A80.png) |
| NLU | NLU.params | `CCCC[C@H](N)C(=O)[O-]` | ![model.name](pngs/NLU.png) |
| Error | 5-hydroxy-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| B96 | beta.beta-diphenyl-alanine.params | `N[C@H](C(=O)[O-])C(c1ccccc1)c1ccccc1` | ![model.name](pngs/B96.png) |
| C27 | homophenylalanine.params | `N[C@@H](CCc1ccccc1)C(=O)[O-]` | ![model.name](pngs/C27.png) |
| B54 | alpha-methyl-3-hydroxy-tyrosine.params | `C[C@H](N)C(=O)[O-].Cc1ccc(O)c(O)c1` | ![model.name](pngs/B54.png) |
| C30 | homoserine.params | `N[C@@H](CCO)C(=O)[O-]` | ![model.name](pngs/C30.png) |
| B61 | alpha-methyl-tyrosine.params | `C[C@H](N)C(=O)[O-].Cc1ccc(O)cc1` | ![model.name](pngs/B61.png) |
| BMAA | 2-amino-3-methylamino-propanoic_acid.params | `C[NH2+]C[C@H](N)C(=O)[O-]` | ![model.name](pngs/BMAA.png) |
| C04 | beta-hydroxy-norvaline.params | `CCC(O)[C@H](N)C(=O)[O-]` | ![model.name](pngs/C04.png) |
| Error | HIP.params | X | (ValueError)  causes seg fault |
| Error | 3-methyl-histidine.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| A30 | 2-amino-4-bromo-4-pentenoic_acid.params | `C=C(Br)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/A30.png) |
| V04 | V04.params | `N[C@H](C(=O)[O-])C(O)c1ccc(O)c(Cl)c1` | ![model.name](pngs/V04.png) |
| Error | 5-fluoro-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| DAB | 2.4-diaminobutyric_acid.params | `N[C@@H](CC[NH3+])C(=O)[O-]` | ![model.name](pngs/DAB.png) |
| C89 | 4-fluoro-proline_puck.params_rot | `CC(F)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/C89.png) |
| HLU | HLU.params | `CC(C)CC[C@H](N)C(=O)[O-]` | ![model.name](pngs/HLU.png) |
| A32 | 2-amino-octanoic_acid.params | `CCCCCC[C@H](N)C(=O)[O-]` | ![model.name](pngs/A32.png) |
| B50 | alpha-amino-glycine.params | `NC([NH3+])C=O` | ![model.name](pngs/B50.png) |
| B95 | beta-beta-dicyclohexyl-alanine.params | `N[C@H](C(=O)[O-])C(C1CCCCC1)C1CCCCC1` | ![model.name](pngs/B95.png) |
| A43 | 2-hydroxy-phenylalanine.params | `N[C@@H](Cc1ccccc1O)C(=O)[O-]` | ![model.name](pngs/A43.png) |
| Error | S56.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| C93 | hexafluoro-leucine.params | `N[C@@H](CC(C(F)(F)F)C(F)(F)F)C(=O)[O-]` | ![model.name](pngs/C93.png) |
| C20 | ethionine.params | `CCSCC[C@H](N)C(=O)[O-]` | ![model.name](pngs/C20.png) |
| B74 | beta-(2-naphthyl)-alanine.params | `N[C@@H](Cc1ccc2ccccc2c1)C(=O)[O-]` | ![model.name](pngs/B74.png) |
| C55 | tert-butyl-glycine.params | `CC(C)(C)[C@H](N)C(=O)[O-]` | ![model.name](pngs/C55.png) |
| B44 | 9-anthryl-alanine.params | `N[C@@H](Cc1c2ccccc2cc2ccccc12)C(=O)[O-]` | ![model.name](pngs/B44.png) |
| Error | 6-chloro-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| C92 | fluoro-leucine_ent2.params | `CC(CF)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/C92.png) |
| B62 | alpha-methyl-valine.params | `CCC.C[C@H](N)C(=O)[O-]` | ![model.name](pngs/B62.png) |
| C41 | penicillamine.params | `CC(C)(S)[C@H](N)C(=O)[O-]` | ![model.name](pngs/C41.png) |
| ABA | ABA.params | `CC[C@H](N)C(=O)[O-]` | ![model.name](pngs/ABA.png) |
| Error | 4-phenyl-phenylalanine_tyr_rot.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| A98 | 4-amino-piperidine-4-carboxylic-acid.params | `CC[NH2+]CC[C@H](N)C(=O)[O-]` | ![model.name](pngs/A98.png) |
| A83 | 3-methyl-histidine_prot.params | `CN1CNC=C1C[C@H](N)C(=O)[O-]` | ![model.name](pngs/A83.png) |
| B59 | alpha-methyl-proline.params | `CCC.C[C@H](N)C(=O)[O-]` | ![model.name](pngs/B59.png) |
| B56 | alpha-methyl-histidine.params | `C[C@H](N)C(=O)[O-].Cc1c[nH]cn1` | ![model.name](pngs/B56.png) |
| A45 | 2-indanyl-glycine_puck2.params | `N[C@H](C(=O)[O-])C1Cc2ccccc2C1` | ![model.name](pngs/A45.png) |
| Error | MPH.params | X | (RuntimeError)  <br/><br/>File: /Volumes/MacintoshHD3/benchmark/W.fujii.release/rosetta.Fujii.release/_commits_/main/source/src/core/conformation/Residue.cc:1365<br/>[ ERROR ] UtilityExitException<br/>ERROR: Unable to fill in missing atoms.<br/><br/> |
| C15 | diphenylglycine.params | `N[C@H](C(=O)[O-])c1ccccc1.c1ccccc1` | ![model.name](pngs/C15.png) |
| Error | 6-methyl-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| Error | 6-bromo-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| BB8 | phenyl-serine.params | `N[C@H](C(=O)[O-])C(O)c1ccccc1` | ![model.name](pngs/BB8.png) |
| Error | 7-azatryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| C54 | tert-butyl-cysteine.params | `CC(C)(C)SC[C@H](N)C(=O)[O-]` | ![model.name](pngs/C54.png) |
| V03 | V03.params | `N[C@H](C(=O)[O-])c1cc(O)cc(O)c1` | ![model.name](pngs/V03.png) |
| MPA | MPA.params | `Cc1ccc(cc1)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/MPA.png) |
| TES | TES.params | `CC(=O)NC1C(NC(=O)C[C@H](N)C(=O)[O-])OC(CO)C(OC2OC(CO)C(OC3OC(COC4OC(COC5OC(CO)C(O)C(O)C5O)C(O)C(OC5OC(CO)C(O)C(O)C5O)C4O)C(O)C(OC4OC(CO)C(O)C(O)C4O)C3O)C(O)C2NC(C)=O)C1O` | ![model.name](pngs/TES.png) |
| C89 | 4-fluoro-proline_puck.params | `CC(F)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/C89.png) |
| C53 | tert-butyl-alanine.params | `CC(C)(C)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/C53.png) |
| A06 | 1-methyl-histidine.params | `Cn1cnc(c1)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/A06.png) |
| A69 | 3-amino-tyrosine.params | `Nc1cc(ccc1O)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/A69.png) |
| HTY | HTY.params | `N[C@@H](Cc1ccc(O)cc1O)C(=O)[O-]` | ![model.name](pngs/HTY.png) |
| B02 | 4-amino-tetrahydropyran-4-carboxylic_acid.params | `CCOCC[C@H](N)C(=O)[O-]` | ![model.name](pngs/B02.png) |
| Error | alpha-methyl-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| HPR | HPR.params | `CC(O)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/HPR.png) |
| Error | 4-methyl-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| C60 | trifluoro-alanine.params | `N[C@H](C(=O)[O-])C(F)(F)F` | ![model.name](pngs/C60.png) |
| C01 | beta-cyclopentyl-alanine.params | `N[C@@H](CC1CCCC1)C(=O)[O-]` | ![model.name](pngs/C01.png) |
| A12 | 2.4-dimethyl-phenylalanine.params | `Cc1ccc(C[C@H](N)C(=O)[O-])c(C)c1` | ![model.name](pngs/A12.png) |
| A07 | 1-methyl-histidine_prot.params | `CN1C=C(C[C@H](N)C(=O)[O-])NC1` | ![model.name](pngs/A07.png) |
| Error | MTP.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| A68 | 3-aminomethyl-phenylalanine.params | `N[C@@H](Cc1cccc(c1)C[NH3+])C(=O)[O-]` | ![model.name](pngs/A68.png) |
| V01 | V01.params | `N[C@H](C(=O)[O-])C(O)c1ccc(O)c(Cl)c1` | ![model.name](pngs/V01.png) |
| C12 | cyclohexyl-glycine.params | `N[C@H](C(=O)[O-])C1CCCCC1` | ![model.name](pngs/C12.png) |
| A84 | 3-methyl-phenylalanine.params | `Cc1cccc(c1)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/A84.png) |
| C61 | trifluoro-leucine.params | `CC(C[C@H](N)C(=O)[O-])C(F)(F)F` | ![model.name](pngs/C61.png) |
| BCS | BCS.params | `N[C@@H](CSCc1ccccc1)C(=O)[O-]` | ![model.name](pngs/BCS.png) |
| A78 | 3-hydroxy-phenylalanine.params | `N[C@@H](Cc1cccc(O)c1)C(=O)[O-]` | ![model.name](pngs/A78.png) |
| C95 | 3-chloro-phenylalanine.params | `N[C@@H](Cc1cccc(Cl)c1)C(=O)[O-]` | ![model.name](pngs/C95.png) |
| A91 | 4.5-dehydro-leucine.params | `C=C(C)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/A91.png) |
| Error | alpha-aminoadipic_acid.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| Error | 7-bromo-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| B63 | amino-ethyl-cysteine.params | `N[C@@H](CSCC[NH3+])C(=O)[O-]` | ![model.name](pngs/B63.png) |
| B19 | 4-fluoro-proline.params | `CC(F)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/B19.png) |
| C03 | beta-fluoro-alanine.params | `N[C@@H](CF)C(=O)[O-]` | ![model.name](pngs/C03.png) |
| C91 | fluoro-leucine_ent1.params | `CC(CF)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/C91.png) |
| A94 | 4-aminomethyl-phenylalanine.params | `N[C@@H](Cc1ccc(cc1)C[NH3+])C(=O)[O-]` | ![model.name](pngs/A94.png) |
| ORN | ornithine.params | `N[C@@H](CCC[NH2+])C(=O)[O-]` | ![model.name](pngs/ORN.png) |
| Error | 5-methyl-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| Error | 5-chloro-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| Error | HP2.params | X | (ValueError)  causes seg fault |
| Error | 4-phenyl-phenylalanine.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| B21 | 4-hydroxy-phenylglycine.params | `N[C@H](C(=O)[O-])c1ccc(O)cc1` | ![model.name](pngs/B21.png) |
| DPP | 2.3-diaminopropionic_acid.params | `N[C@@H](C[NH3+])C(=O)[O-]` | ![model.name](pngs/DPP.png) |
| IGL | 2-indanyl-glycine_puck1.params | `N[C@H](C(=O)[O-])C1Cc2ccccc2C1` | ![model.name](pngs/IGL.png) |
| C05 | beta-iodo-alanine.params | `N[C@@H](CI)C(=O)[O-]` | ![model.name](pngs/C05.png) |
| A04 | 1-amino-cyclopentane-carboxylic_acid.params | `CCCC[C@H](N)C(=O)[O-]` | ![model.name](pngs/A04.png) |
| A92 | 4.5-dehydro-lysine.params | `N[C@@H](CC=CC[NH3+])C(=O)[O-]` | ![model.name](pngs/A92.png) |
| A48 | 2-methyl-phenylalanine.params | `Cc1ccccc1C[C@H](N)C(=O)[O-]` | ![model.name](pngs/A48.png) |
| Error | SAL.params | X | (RuntimeError)  <br/><br/>File: /Volumes/MacintoshHD3/benchmark/W.fujii.release/rosetta.Fujii.release/_commits_/main/source/src/core/conformation/Residue.cc:200<br/>[ ERROR ] UtilityExitException<br/>ERROR: <br/><br/> |
| B67 | beta-(1-naphthyl)-alanine.params | `N[C@@H](Cc1cccc2ccccc12)C(=O)[O-]` | ![model.name](pngs/B67.png) |
| C42 | phenylglycine.params | `N[C@H](C(=O)[O-])c1ccccc1` | ![model.name](pngs/C42.png) |
| B31 | 4-tert-butyl-phenylalanine.params | `CC(C)(C)c1ccc(cc1)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/B31.png) |
| APA | APA.params | `Nc1ccc(cc1)C[C@H](N)C(=O)[O-]` | ![model.name](pngs/APA.png) |
| Error | 6-fluoro-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| Error | 4-carboxy-phenylalanine.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| Error | BZP.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| Error | n-in-methyl-tryptophan.params | X | (ArgumentError)  Python argument types in<br/>    rdkit.Chem.rdmolops.ReplaceSubstructs(NoneType, Mol, Mol)<br/>did not match C++ signature:<br/>    ReplaceSubstructs(RDKit::ROMol mol, RDKit::ROMol query, RDKit::ROMol replacement, bool replaceAll=False, unsigned int replacementConnectionPoint=0, bool useChirality=False) |
| C16 | dipropyl-glycine.params | `CCC.CCC[C@H](N)C(=O)[O-]` | ![model.name](pngs/C16.png) |