import os
import shutil
from ase import Atoms
from ase.units import kcal, mol
from ase.io import read
from fennol.ase import FENNIXCalculator

def FENNIX_OPT(poltype,inputstruct,resfilename):
    cmdstr = f"geometric-optimize --ase-class=fennol.ase.FENNIXCalculator --ase-kwargs=$CONFIG --engine=ase {inputstruct} {resfilename}"
    poltype.call_subsystem([cmdstr], True)
    opted_xyz = inputstruct.replace('.xyz', '_optim.xyz') 
    if os.path.isfile(opted_xyz):
      lines = open(opted_xyz).readlines()
      natom = int(lines[0].split()[0])

    prefix = resfilename.split('_constr.txt')[0]
    fname = inputstruct.split(".xyz")[0]
    shutil.move(f"{fname}.log", f"{prefix}.log")
    optxyz = f"{prefix}.xyz"
    with open(optxyz, 'w') as f:
      for i in range(len(lines) - (natom+2), len(lines)):
        f.write(lines[i])

    return optxyz, f"echo 'Doing FENNIX optimization for {prefix}' "

def FENNIX_SP(poltype,optxyz,output):
  atoms = read(optxyz)
  atom_symbols_list = atoms.get_chemical_symbols()
  positions_2Darray = atoms.positions
  
  # positions have to be in Angstrom
  molecule = Atoms(
      symbols=atom_symbols_list, 
      positions=positions_2Darray
      )
  
  #calc = FENNIXCalculator(model='$FENNIX_MODEL_PATH')
  calc = FENNIXCalculator(model='/home/liuchw/TryOtherSoftware/FENNIX/ani2x_model0.fnx')
  molecule.calc = calc
  
  # ASE returns energies in eV (can be converted using units from ase.units module)
  fennix_energy = molecule.get_potential_energy() 
  
  # convert from eV to kcal/mol (it's inverted intentionally)
  fennix_energy = fennix_energy * mol / kcal
  
  with open(output, 'w') as f:
    f.write(f'Normal Termination Calling FENNIX_SP: {optxyz}\n')
    f.write(f'Total Energy: {fennix_energy:.4f} Kcal/mol\n')
  return 