#
# This program is used to find the diheral angles that are eligible 
# for restrained optimization. The results can be used for torsion restrained
# optimzation both in _opt-1. file or torsion -opt file
# 

# Author: Chengwen Liu
# Date: Jul. 2024

from rdkit import Chem

def findTorsionToBeRestrained(mol):
  nh2_smt = "[NH2]([H])([H])[*]"
  pattern = Chem.MolFromSmarts(nh2_smt)
  matches = mol.GetSubstructMatches(pattern)
  HsOfNitrogen = {}
  for match in matches:
    n, h1, h2, _ = list(match)
    HsOfNitrogen[n] = [h1, h2]
  
  rotbond2torsion = {}
  for bond in mol.GetBonds():
    tmp = []
    bond_type = bond.GetBondType()
    # only SINGLE bond is considered rotatable
    if (bond_type == Chem.rdchem.BondType.SINGLE):
      b1 = bond.GetBeginAtomIdx()
      b2 = bond.GetEndAtomIdx()
      left_index = -1 
      right_index = -1 
      left_mass = 0 
      right_mass = 0 
      
      left_trigonal = []
      right_trigonal = []
      
      if b1 in HsOfNitrogen.keys():
        left_trigonal = HsOfNitrogen[b1]
      if b2 in HsOfNitrogen.keys():
        right_trigonal = HsOfNitrogen[b2]

      # get the heaviest atom of left
      for neigh in mol.GetAtomWithIdx(b1).GetNeighbors():
        neig_idx = neigh.GetIdx()
        if neig_idx != b2:
          atomicNum = mol.GetAtomWithIdx(neig_idx).GetAtomicNum()
          if atomicNum > left_mass:
            left_mass = atomicNum
            left_index = neig_idx
          
          if (atomicNum == left_mass) and (neig_idx < left_index):
            left_index = neig_idx
              
      # get the heaviest atom of right 
      for neigh in mol.GetAtomWithIdx(b2).GetNeighbors():
        neig_idx = neigh.GetIdx()
        if neig_idx != b1:
          atomicNum = mol.GetAtomWithIdx(neig_idx).GetAtomicNum()
          if atomicNum > right_mass:
            right_mass = atomicNum
            right_index = neig_idx
          
          if (atomicNum == right_mass) and (neig_idx < right_index):
            right_index = neig_idx
      
      if (left_trigonal != []) and (right_index != -1):
        for l in left_trigonal:
          if ([l, b1, b2, right_index] not in tmp):
            if f'{b1} {b2}' not in rotbond2torsion.keys():
              rotbond2torsion[f'{b1} {b2}'] = [[l, b1, b2, right_index]]
            else:
              rotbond2torsion[f'{b1} {b2}'] += [[l, b1, b2, right_index]]
              
            tmp.append([l, b1, b2, right_index])
      
      if (right_trigonal != []) and (left_index != -1):
        for r in right_trigonal: 
          if ([left_index, b1, b2, r] not in tmp):
            if f'{b1} {b2}' not in rotbond2torsion.keys():
              rotbond2torsion[f'{b1} {b2}'] = [[left_index, b1, b2, r]]
            else:
              rotbond2torsion[f'{b1} {b2}'] += [[left_index, b1, b2, r]]
            tmp.append([left_index, b1, b2, r])
      
      if (left_index != -1) and (right_index != -1) and ([left_index, b1, b2, right_index] not in tmp):
        rotbond2torsion[f'{b1} {b2}'] = [[left_index, b1, b2, right_index]]
        tmp.append([left_index, b1, b2, right_index])
   
  return rotbond2torsion
