import shutil
from rdkit import Chem

# this is to modify the intermediate tinker key files generated by poltype

def modkey2(poltype):
  rdkitmol = Chem.MolFromMolFile(poltype.molstructfname,removeHs=False)
  scaleDipoleAtoms = []
  fixDipoleAtoms = []
  # acid or ester
  smt1 = "[C](=O)[O][#6,#1]" 
  # acetate anion 
  smt2 = "[C](=O)[O-]" 
  smts = [smt1, smt2]
  for smt in smts:
    pattern = Chem.MolFromSmarts(smt)
    match = rdkitmol.GetSubstructMatches(pattern)
    if match:
      for i in range(len(match)):	
        scaleDipoleAtoms += match[i][1:3]
        fixDipoleAtoms += match[i][0:3]
  
  # scale dipole of the atoms
  scaleDipoleRatio = 0.3
  atom2type = {}
  lines = open(poltype.xyzoutfile).readlines()
  for line in lines[1:]:
    ss = line.split()
    atom2type[int(ss[0])] = ss[5]

  scaleDipoleTypes = [atom2type[int(i)+1] for i in scaleDipoleAtoms]
  
  tmpkeyfile = poltype.key2fnamefromavg + '_tmp'
  lines = open(poltype.key2fnamefromavg).readlines()
  lines_append = []
  for atm in fixDipoleAtoms:
    lines_append.append(f"FIX-ATOM-DIPOLE {int(atm)+1} X\n")
    lines_append.append(f"FIX-ATOM-DIPOLE {int(atm)+1} Y\n")
    lines_append.append(f"FIX-ATOM-DIPOLE {int(atm)+1} Z\n")
  
  for i in range(len(lines)):
    line = lines[i]
    if ('multipole ' in line) and (line.split()[1] in scaleDipoleTypes):
      [dx, dy, dz] = lines[i+1].split()
      lines[i+1] = ' '*37 + f"{float(dx)*scaleDipoleRatio:10.5f}{float(dy)*scaleDipoleRatio:10.5f}{float(dz)*scaleDipoleRatio:10.5f}\n"
      
  with open(tmpkeyfile, 'w') as f:
    for line in lines:
      if "RESP-WEIGHT " in line:
        f.write(line)
        for apline in lines_append:
          f.write(apline)
          f.write(apline)
          f.write(apline)
      else:
        f.write(line)

  # rename key2 to key2b
  shutil.move(poltype.key2fnamefromavg, poltype.key2fnamefromavg+'b')
  # rename key2_tmp to key2
  shutil.move(poltype.key2fnamefromavg + '_tmp', poltype.key2fnamefromavg)