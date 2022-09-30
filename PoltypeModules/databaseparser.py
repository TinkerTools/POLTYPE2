from tqdm import tqdm
import warnings
import math
import sys
import os
from packaging import version
import rdkit
from rdkit.Chem import rdFMCS
from openbabel import openbabel
from rdkit.Chem import rdmolfiles
import itertools
import re
from rdkit import Chem
import copy
import symmetry as symm
from rdkit.Chem.Lipinski import RotatableBondSmarts
from rdkit.Chem.rdMolDescriptors import CalcNumRings
import numpy as np
import json
import torsionfit as torfit
from rdkit import DataStructs
import torsiongenerator as torgen
from itertools import combinations
import shutil
 

def CheckIfStringIsFloat(string):
    """
    Intent: Check if string is float value
    Input: String
    Output: boolean indicating if string is a float or not
    Referenced By: appendtofile
    Description: 
    1. Assume is not float
    2. Try to convert to float, if works set boolean to True
    3. If fail return boolean
    """
    # STEP 1
    isfloat=False
    # STEP 2
    try:
        float(string)
        isfloat=True
    except:
        # STEP 3
        pass
    return isfloat


def appendtofile(poltype, vf,newname, bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo,polarprmstotransferinfo,soluteprms,amoebaplusvdwprmstotransferinfo,ctprmstotransferinfo,cpprmstotransferinfo,bondcfprmstotransferinfo,anglecfprmstotransferinfo,tortorprmstotransferinfo):
    """
    Intent: Append the parameters found from database search to key file
    Input: Original key file, new name of updated key file, various dictionaries of parameters to comments about 
    Output: -  
    Referenced By: poltype.py - GenerateParameters 
    Description: 
    1. Iterate over lines of key files 
    2. Detect when find atom block and when after atom block
    3. If polarize or multipole in line, skip this line (only adding valence and vdw parameters for AMOEBA from database here).
    4. If AMOEBA+ forcefield, write out AMOEBA+ vdw, else write out AMOEBA vdw
    5. Write out bond parameters
    6. Write out angle parameters
    7. Write out stretch bend parameters
    8. Write out op-bend parameters
    9. Write out torsion parameters
    10.Write out solute parameters
    11.Write out tor-tor parameters
    12.If AMOEBA+ forcefield, write out charge penetration, charge transfer and charge flux parameters
    13.If user gives key file parameters as input to write in, then write those to key file
    """
    tempname=vf.replace('.key','_temp.key')
    f=open(tempname,'w')
    if poltype.writeoutpolarize==True and poltype.writeoutmultipole==True: 
        temp=open(vf,'r')
        results=temp.readlines()
        temp.close()
        foundatomblock=False
        atomline=False
        wroteout=False
        linestoskip=[]
        # STEP 1
        for theline in results:
            linesplit=theline.split()
            # STEP 2
            if 'atom' in theline:
                atomline=True
                if foundatomblock==False:
                    foundatomblock=True

                        
            else:
                atomline=False
            # STEP 3
            if 'polarize' in theline and poltype.writeoutpolarize==False:
                linestoskip.append(theline)
            if 'multipole' in theline and poltype.writeoutmultipole==False:
                linestoskip.append(theline)
            if len(linesplit)>0:
                if CheckIfStringIsFloat(linesplit[0])==True and poltype.writeoutmultipole==False:
                    linestoskip.append(theline)


            if foundatomblock==True and atomline==False and wroteout==False:
                wroteout=True
                f.write('\n')
                # STEP 4
                if poltype.writeoutvdw==True:
                    if poltype.forcefield=='AMOEBA+':
                        for line,transferinfo in amoebaplusvdwprmstotransferinfo.items():
                            f.write(transferinfo)
                            f.write(line)
                            f.write('\n')

                    else:
                        for line,transferinfo in vdwprmstotransferinfo.items():
                            f.write(transferinfo)
                            f.write(line)
                            f.write('\n')
                
                # STEP 5
                if poltype.writeoutbond==True:
                    for line,transferinfo in bondprmstotransferinfo.items():
                        f.write(transferinfo)
                        f.write(line)
                        f.write('\n')
                # STEP 6 
                if poltype.writeoutangle==True:
                    for line,transferinfo in angleprmstotransferinfo.items():
                        f.write(transferinfo)
                        f.write(line)
                        f.write('\n')
                # STEP 7
                if poltype.writeoutstrbnd==True:
                    for line,transferinfo in strbndprmstotransferinfo.items():
                        f.write(transferinfo)
                        f.write(line)
                        f.write('\n')
                # STEP 8
                if poltype.writeoutopbend==True:
                    for line,transferinfo in opbendprmstotransferinfo.items():
                        f.write(transferinfo)
                        f.write(line)
                        f.write('\n')
                # STEP 9
                if poltype.writeouttorsion==True:
                    for line,transferinfo in torsionprmstotransferinfo.items():
                        f.write(transferinfo)
                        f.write(line)
                        f.write('\n')
                # STEP 10
                for line in soluteprms:
                    f.write(line)
                    f.write('\n')
                # STEP 11
                for line,transferinfo in tortorprmstotransferinfo.items():
                    if 'tortors' in line:
                        f.write(transferinfo)
                    f.write(line)
                    f.write('\n')
                # STEP 12
                if poltype.forcefield=='AMOEBA+':
                    for line,transferinfo in ctprmstotransferinfo.items():
                        f.write(transferinfo)
                        f.write(line)
                        f.write('\n')
                    for line,transferinfo in cpprmstotransferinfo.items():
                        f.write(transferinfo)
                        f.write(line)
                        f.write('\n')
                    for line,transferinfo in bondcfprmstotransferinfo.items():
                        f.write(transferinfo)
                        f.write(line)
                    f.write('\n')
                    for line,transferinfo in anglecfprmstotransferinfo.items():
                        f.write(transferinfo)
                        f.write(line)
                        f.write('\n')


                                        
            else:
                if theline not in linestoskip:
                    f.write(theline)
    
    # STEP 13
    if poltype.inputkeyfile!=None:
        if poltype.writeouttorsion==True:
            for line,transferinfo in torsionprmstotransferinfo.items():
                f.write(transferinfo)
                f.write(line)
                f.write('\n')
        temp=open(poltype.inputkeyfile,'r')
        results=temp.readlines()
        temp.close()
        f.write("parameters " + poltype.paramhead + "\n")
        for line in results: # handle case where user gives giant prm file as input
            linesplit=line.split()
            if len(linesplit)>1:
                if linesplit[0]=='torsion':
                    if poltype.writeouttorsion==True:
                        for theline,transferinfo in torsionprmstotransferinfo.items():
                            thelinesplit=theline.split()
                            tors=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])]
                            thetors=[int(thelinesplit[1]),int(thelinesplit[2]),int(thelinesplit[3]),int(thelinesplit[4])]
                            if tors==thetors or tors==thetors[::-1]:
                                line=theline
                                break
                    
            f.write(line)

    f.close()
    os.rename(tempname,newname)

def ReadSmallMoleculeLib(poltype,filepath):
    """
    Intent: Read amoeba09 SMARTS, atom order -> amoeba09 tinker description file
    Input: Filepath
    Output: Dictionary of SMARTS, atom order -> amoeba09 tinker description
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over lines of file  
    2. Grab symbol
    3. Grab tinker description
    4. Grab SMARTS string
    5. Grab atom order list
    6. Put items in dictionary 
    """
    temp=open(filepath,'r')
    results=temp.readlines()
    temp.close()
    smartsatomordertoelementtinkerdescrip={}
    # STEP 1
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if len(linesplit)==0:
            continue
        # STEP 2
        elementsymb=linesplit[0]
        newline=' '.join(linesplit[1:])
        newsplit=newline.split('%')
        # STEP 3
        tinkerdescrip=newsplit[0].lstrip().rstrip()
        # STEP 4
        smarts=newsplit[1].lstrip().rstrip()
        # STEP 5
        atomindices=newsplit[2].lstrip().rstrip()
        atomorderlist=atomindices.split()
        atomorderlist=tuple([int(i) for i in atomorderlist])
        ls=[smarts,atomorderlist]
        newls=[elementsymb,tinkerdescrip]
        # STEP 6
        smartsatomordertoelementtinkerdescrip[tuple(ls)]=tuple(newls)
    return smartsatomordertoelementtinkerdescrip

def GrabParameters(poltype,fname):
    """
    Intent: Grab parameters from parameter file
    Input:
    Output:
    Referenced By: mutation.py GrabBgnToEndPrms
    Description: 
    1. Iterate over lines of key file
    2. If atom keyword, save atom line
    3. If bond keyword, save bond parameter line
    4. If angle keyword, save angle parameter line
    5. If torsion keyword, save torsion parameter line
    6. If strbnd keyword, save stretch bend parameter line
    7. If opbend keyword, save opbend parameter line
    8. If polarize keyword, save polarize parameter line
    9. If vdw keyword, save vdw parameter line
    10. If multipole keyword, save multipole parameter lines 
    """
    atomdefs=[]
    bondprms=[]
    angleprms=[]
    torsionprms=[]
    strbndprms=[]
    opbendprms=[]
    polarizeprms=[]
    vdwprms=[]
    mpoleprms=[]
    temp=open(fname,'r')
    results=temp.readlines()
    temp.close()
    # STEP 1
    for lineidx in range(len(results)):
        line=results[lineidx]
        # STEP 2
        if 'atom' in line:
            atomdefs.append(line)
        # STEP 3
        elif 'bond' in line:
            bondprms.append(line)
        # STEP 4
        elif 'angle' in line or 'anglep' in line:
            angleprms.append(line)
        # STEP 5
        elif 'torsion' in line:
            torsionprms.append(line) 
        # STEP 6
        elif 'strbnd' in line:
            strbndprms.append(line)
        # STEP 7
        elif 'opbend' in line:
            opbendprms.append(line)
        # STEP 8
        elif 'polarize' in line:
            polarizeprms.append(line)
        # STEP 9
        elif 'vdw' in line: 
            vdwprms.append(line)
        # STEP 10
        elif 'multipole' in line:
            mpolelist=[line,results[lineidx+1],results[lineidx+2],results[lineidx+3],results[lineidx+4]]
            for mpoleline in mpolelist:
                mpoleprms.append(mpoleline)

    return atomdefs,bondprms,angleprms,torsionprms,strbndprms,opbendprms,polarizeprms,vdwprms,mpoleprms
 
def GrabParametersFromPrmFile(poltype,bondtinkerclassestopoltypeclasses,opbendtinkerclassestopoltypeclasses,opbendtinkerclassestotrigonalcenterbools,angletinkerclassestopoltypeclasses,torsiontinkerclassestopoltypeclasses,poltypetoprmtype,atomtinkerclasstopoltypeclass,typestoframedefforprmfile,fname,skipmultipole=False):
    """
    Intent: Grab parameters from amoeba09 parameter file.
    Input: dictionaries of tinker class -> poltype class numbers for each parameter type
    Output: arrays of parameter lines for each parameter type in amoeba09
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over lines of parameter file
    2. If bond keyword in line, check if bond classes are in input dictionaries and if so, save parameter line 
    3. If angle keyword in line, check if angle classes are in input dictionaries and if so, save parameter line 
    4. If torsion keyword in line, check if torsion classes are in input dictionaries and if so, save parameter line 
    5. If strbnd keyword in line, check if strbnd classes are in input dictionaries and if so, save parameter line 
    6. If opbend keyword in line, check if opbend classes are in input dictionaries and if so, save parameter line 
    7. If multipole keyword in line, check if multipole classes are in input dictionaries and if so, save parameter line 
    8. If polarize keyword in line, check if polarize classes are in input dictionaries and if so, save parameter line 
    9. If vdw keyword in line, check if vdw classes are in input dictionaries and if so, save parameter line 
    10. If pitor keyword in line, check if middle two classes are in torsion input dictionaries, if so map parameter line to torsion parameters (later on will modify torsion parameters to account for not adding pitor parameters to keyfile)
    """
    temp=open(fname,'r') 
    results=temp.readlines()
    temp.close()
    bondprms=[]
    angleprms=[]
    torsionprms=[]
    strbndprms=[]
    mpoleprms=[]
    opbendprms=[]
    polarizeprms=[]
    vdwprms=[]
    torsiontopitor={}
    # STEP 1
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        linesplitall=re.split(r'(\s+)', line)
        if '#' in line:
            continue
        # STEP 2
        if 'bond' in line and 'cubic' not in line and 'quartic' not in line:
            bondclasslist=[int(linesplit[1]),int(linesplit[2])]
            foundbond=False
            if tuple(bondclasslist) in bondtinkerclassestopoltypeclasses.keys():
                bondtup=tuple(bondclasslist)
                foundbond=True
            elif tuple(bondclasslist[::-1]) in bondtinkerclassestopoltypeclasses.keys():
                bondtup=tuple(bondclasslist[::-1])
                foundbond=True
            if foundbond==True:
                classes=bondtinkerclassestopoltypeclasses[bondtup]
                for boundcls in classes:
                    linesplitall[2]=str(boundcls[0])    
                    linesplitall[4]=str(boundcls[1])  
                    newline=''.join(linesplitall)
                    bondprms.append(newline)    
        # STEP 3       
        elif 'angle-cubic' not in line and 'angle-quartic' not in line and 'pentic' not in line and 'sextic' not in line and ('angle' in line or 'anglep' in line) :
            angleclasslist=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
            foundangle=False
            if tuple(angleclasslist) in angletinkerclassestopoltypeclasses.keys():
                angletup=tuple(angleclasslist)
                foundangle=True
            elif tuple(angleclasslist[::-1]) in angletinkerclassestopoltypeclasses.keys():
                angletup=tuple(angleclasslist[::-1])
                foundangle=True
            if foundangle==True:
                classes=angletinkerclassestopoltypeclasses[angletup]
                for boundcls in classes:
                    linesplitall[2]=str(boundcls[0])    
                    linesplitall[4]=str(boundcls[1])  
                    linesplitall[6]=str(boundcls[2])
                    newline=''.join(linesplitall)
                    angleprms.append(newline) 
        # STEP 4 
        elif 'torsion' in line and 'torsionunit' not in line:
            torsionclasslist=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])]
            foundtorsion=False
            if tuple(torsionclasslist) in torsiontinkerclassestopoltypeclasses.keys():
                torsiontup=tuple(torsionclasslist)
                foundtorsion=True
            elif tuple(torsionclasslist[::-1]) in torsiontinkerclassestopoltypeclasses.keys():
                torsiontup=tuple(torsionclasslist[::-1])
                foundtorsion=True
            if foundtorsion==True:   
                classes=torsiontinkerclassestopoltypeclasses[torsiontup]
                for boundcls in classes:
                    if linesplitall[0]=='': # sometimes when torsion is added back to the key file, it has a space in front of it
                        linesplitall[4]=str(boundcls[0])    
                        linesplitall[6]=str(boundcls[1])  
                        linesplitall[8]=str(boundcls[2])
                        linesplitall[10]=str(boundcls[3])
                    else:
                        linesplitall[2]=str(boundcls[0])    
                        linesplitall[4]=str(boundcls[1])  
                        linesplitall[6]=str(boundcls[2])
                        linesplitall[8]=str(boundcls[3])

                    newline=''.join(linesplitall)
                    torsionprms.append(newline)
        # STEP 5
        elif 'strbnd' in line:
            angleclasslist=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
            foundstrbnd=False
            rev=False
            prm1=linesplitall[8]
            prm2=linesplitall[10]
            if tuple(angleclasslist) in angletinkerclassestopoltypeclasses.keys():
                angletup=tuple(angleclasslist)
                foundstrbnd=True
            elif tuple(angleclasslist[::-1]) in angletinkerclassestopoltypeclasses.keys():
                angletup=tuple(angleclasslist[::-1])
                foundstrbnd=True
                rev=True
            if foundstrbnd==True:
                classes=angletinkerclassestopoltypeclasses[angletup]
                for boundcls in classes:
                    linesplitall[2]=str(boundcls[0])    
                    linesplitall[4]=str(boundcls[1])  
                    linesplitall[6]=str(boundcls[2])
                    if rev==True:
                        linesplitall[8]=prm2
                        linesplitall[10]=prm1
                    newline=''.join(linesplitall)
                    strbndprms.append(newline) 

               
        # STEP 6 
        elif 'opbend' in line and 'opbendtype' not in line and 'cubic' not in line and 'quartic' not in line and 'pentic' not in line and 'sextic' not in line:
            bondclasslist=[int(linesplit[1]),int(linesplit[2])]
            foundopbend=False

            if tuple(bondclasslist) in opbendtinkerclassestopoltypeclasses.keys():
                bondtup=tuple(bondclasslist)
                foundopbend=True

            if foundopbend==True:
                classes=opbendtinkerclassestopoltypeclasses[bondtup]
                boolarray=opbendtinkerclassestotrigonalcenterbools[bondtup]

                for boundcls in classes:
                    linesplitall[2]=str(boundcls[0])    
                    linesplitall[4]=str(boundcls[1])  
                    newlinesplitall=linesplitall[:]
                    newlinesplitall[4]=str(boundcls[0])    
                    newlinesplitall[2]=str(boundcls[1])
                    newline=''.join(linesplitall)
                    if boolarray[1]==True: 
                        opbendprms.append(newline)
        # STEP 7
        elif 'multipole' in line and skipmultipole==False:
            newlinesplit=linesplit[1:-1]
            frames=[int(i) for i in newlinesplit]
            grabit=False
            if tuple(frames) in typestoframedefforprmfile.keys():
                theframe=tuple(frames)
                grabit=True
            elif tuple(frames[::-1]) in typestoframedefforprmfile.keys():
                theframe=tuple(frames[::-1])
                grabit=True
            if grabit==True:
                chgpartofline=linesplitall[-2:] # include space
                phrasepartofline=linesplitall[:2] # include space
                spacetoadd='   '
                newlist=[]
                newlist.extend(phrasepartofline)
                framedef=typestoframedefforprmfile[atomtype]
                for typenum in framedef:
                    newlist.append(str(typenum))
                    newlist.append(spacetoadd)
                newlist=newlist[:-1] # remove last space then add space before chg
                newlist.extend(chgpartofline)
                newline=''.join(newlist)
                mpolelist=[newline,results[lineidx+1],results[lineidx+2],results[lineidx+3],results[lineidx+4]]
                for mpoleline in mpolelist:
                    mpoleprms.append(mpoleline)
        # STEP 8 
        elif 'polarize' in line: 
            atomtype=int(linesplit[1])
            if atomtype in poltypetoprmtype.keys():
                prmtype=poltypetoprmtype[atomtype]
                newline=line.replace('\n','')+' '+str(prmtype)+'\n'
                polarizeprms.append(newline)
        # STEP 9 
        elif 'vdw' in line and 'type' not in line and 'scale' not in line: 
            atomclass=int(linesplit[1])
            atomclasslist=tuple([atomclass])
            if atomclasslist in atomtinkerclasstopoltypeclass.keys():
                prmclasses=atomtinkerclasstopoltypeclass[atomclasslist]
                for prmclasslist in prmclasses:
                    for prmclass in prmclasslist:
                        linesplit[1]=str(prmclass)
                        newline=' '.join(linesplit)+'\n'
                        vdwprms.append(newline)
        # STEP 10
        elif 'pitor' in line:
            linesplit=line.split()
            b=linesplit[1]
            c=linesplit[2]
            for tinkerclasses,poltypeclasses in torsiontinkerclassestopoltypeclasses.items():
                tb=str(tinkerclasses[1])
                tc=str(tinkerclasses[2])
                if (b==tb and c==tc) or (b==tc and c==tb):
                    for cls in poltypeclasses:
                        torsiontopitor[tuple(cls)]=line

    
    return bondprms,angleprms,torsionprms,strbndprms,mpoleprms,opbendprms,polarizeprms,vdwprms,torsiontopitor


def GrabTypeAndClassNumbers(poltype,prmfile):
    """
    Intent: Grab dictionaries of tinker element description -> tinker type and tinker type -> tinker class from amoeba09 parameter file. 
    Input: amoeba09 parameterfile
    Output: dictionaries of tinker element description -> tinker type and tinker type -> tinker class 
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over lines of parameter file, if atom keyword is in line
    2. Grab tinker type
    3. Grab class type
    4. Grab element
    5. Grab tinker description
    6. Put everything into dictionaries 
    """
    temp=open(prmfile,'r')
    results=temp.readlines()
    temp.close()
    elementtinkerdescriptotinkertype={}
    tinkertypetoclass={}
    # STEP 1
    for line in results:
        if 'atom' in line:
            linesplit=line.split()    
            newlinesplit=linesplit[1:-3]
            # STEP 2
            tinkertype=newlinesplit[0]
            # STEP 3
            classtype=newlinesplit[1]
            # STEP 4
            element=newlinesplit[2]
            # STEP 5
            tinkerdescrip=' '.join(newlinesplit[3:])
            ls=[element,tinkerdescrip]
            # STEP 6
            elementtinkerdescriptotinkertype[tuple(ls)]=tinkertype
            tinkertypetoclass[tinkertype]=classtype
    return elementtinkerdescriptotinkertype,tinkertypetoclass



def GrabAtomsForParameters(poltype,mol):
    """
    Intent: Need a list of atoms, bonds, angles and torsions for input molecule to process and parse database later.
    Input: Openbabel mol object
    Output: list of atoms, bonds, angles and torsions
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description: 
    1. Iterate over atoms append atom indices to atom array
    2. Iterate over neighbor of atoms and if bond (atom,natom) hasnt been found append to array
    3. Iterate over neighbor of neighbor and if first atom and neighbor of neighbor are not same atoms, append angle to array
    4. Iterate over neighbors again to find torsion, make sure no indices are repeating, and append torsion to array
    """
    # now we can define arrays to collect bonds, angles and torsions
    listoftorsionsforprm=[]
    listofbondsforprm=[]
    listofanglesforprm=[]
    listofatomsforprm=[]
    # STEP 1
    for atom in openbabel.OBMolAtomIter(mol):
        atomidx=atom.GetIdx()-1
        listofatomsforprm.append([atomidx])
        neighbs=[natom for natom in openbabel.OBAtomAtomIter(atom)]
        # STEP 2
        for natom in neighbs:
            nidx=natom.GetIdx()-1
            bondset=[nidx,atomidx]
            if bondset not in listofbondsforprm and bondset[::-1] not in listofbondsforprm:
                listofbondsforprm.append(bondset)
            # STEP 3
            nextneighbs=[nextatom for nextatom in openbabel.OBAtomAtomIter(natom)]
            for nextneighb in nextneighbs:
                nextneighbidx=nextneighb.GetIdx()-1
                if nextneighbidx!=atomidx:
                    angleset=[nextneighbidx,nidx,atomidx]
                    if angleset not in listofanglesforprm and angleset[::-1] not in listofanglesforprm:
                        listofanglesforprm.append(angleset)
                    nextnextneighbs=[nextnextatom for nextnextatom in openbabel.OBAtomAtomIter(nextneighb)]
                    # STEP 4
                    for nextnextatom in nextnextneighbs:
                        nextnextatomidx=nextnextatom.GetIdx()-1
                        if nextnextatomidx!=nidx and nextnextatomidx!=atomidx:                
                            torsionset=[nextnextatomidx,nextneighbidx,nidx,atomidx]
                            if torsionset not in listoftorsionsforprm and torsionset[::-1] not in listoftorsionsforprm:
                                listoftorsionsforprm.append(torsionset)
    return listofatomsforprm,listofbondsforprm,listofanglesforprm,listoftorsionsforprm



def ExtendByNeighbors(poltype,ls):    
    """
    Intent: Add neighbors to input array, when wanting to check if neighbors exist in SMARTS match
    Input: Array of atom indices
    Output: Array of atom indices plus neighbors
    Referenced By: MatchAllPossibleSMARTSToParameterSMARTS
    Description: 
    1. If input array is of length 1 (vdw atom matches) want neighbor of neighbors also (so match neighbor of neighbor on hydrogen for example to get enough information if match is correct)
    2. Iterate over atom indices in input array
    3. Iterate over neighboring atom indices of atom index and append to new array
    4. If input array is of length 1 also iterate over neighbor of neighbors and append indices to new array 
    """
    # STEP 1
    extendneighbofneighb=False
    if len(ls)==1:
        extendneighbofneighb=True 
    newls=[]
    # STEP 2
    for atomidx in ls:
        if atomidx not in newls:
            newls.append(atomidx)
        atom=poltype.rdkitmol.GetAtomWithIdx(atomidx)
        # STEP 3
        for natom in atom.GetNeighbors():
            natomidx=natom.GetIdx()
            if natomidx not in newls:
                newls.append(natomidx)     
            # STEP 4 
            if extendneighbofneighb==True and len(atom.GetNeighbors())==1:
                for nnatom in natom.GetNeighbors():
                    nnatomidx=nnatom.GetIdx()
                    if nnatomidx not in newls:
                        newls.append(nnatomidx)


    return newls
             


def GenerateFragmentSMARTS(poltype,ls):
    """
    Intent: When matching input molecule to database need a way to take part of input molecule, convert to SMARTS to attempt SMARTS matching. This functions generate SMARTS from input array of atom indices in molecule. 
    Input: Array of atom indices from input molecule
    Output: SMART string, second SMARTS that is modified so that any atoms that were aromatic in original mol are aromatic in SMARTS
    Referenced By: GenerateFragmentSMARTSList
    Description:
    1. Create new mol object
    2. Iterate over input array of atom indices
    3. Grab the atom object from input molecule and add atom to newly created mol object
    4. Keep track of old atom index and new atom index in dicttoinary
    5. Keep track if atom is aromatic atom or not in original molecule
    6. Iterate over bonds
    7. If one atom is in new molecule and but the other atom is not, then keep track of this bond as a cut bond. If both are in new molecule then add bond to new molecule.
    8. Generate SMARTS from new mol  
    9. Iterate back over atoms in new mol and if they were aromatic atoms in original molecule, make them aromatic in new mol
    10.Iterate back over bonds in new mol and if they were aromatic bonds in original molecule, make them aromatic in new mol 
    11.Generate new SMARTS from modified mol (containing aromaticity information)
    """
    # STEP 1
    newmol = Chem.Mol()
    mw = Chem.RWMol(newmol)
    # need to treat ring bonds as aromatic since all transferred parameters from amoeba09 are aromatic rings
    oldindextonewindex={}
    aromaticindices=[]
    # STEP 2
    for i,idx in enumerate(ls):
        # STEP 3
        oldatom=poltype.rdkitmol.GetAtomWithIdx(idx)
        mw.AddAtom(oldatom)
        # STEP 4
        oldindextonewindex[idx]=i
        oldatombabel=poltype.mol.GetAtom(idx+1)
        # STEP 5
        isaro=oldatombabel.IsAromatic()
        isinring=oldatombabel.IsInRing()
        hyb=oldatombabel.GetHyb()
        if isaro==True and isinring==True and hyb==2:
            aromaticindices.append(i)

    atomswithcutbonds=[]
    aromaticbonds=[]
    bonditer=poltype.rdkitmol.GetBonds()
    # STEP 6
    for bond in bonditer:
        oendidx = bond.GetEndAtomIdx()
        obgnidx = bond.GetBeginAtomIdx()
        babelbond=poltype.mol.GetBond(oendidx+1,obgnidx+1)
        isinring=babelbond.IsInRing()
        # STEP 7
        if oendidx in oldindextonewindex.keys() and obgnidx not in oldindextonewindex.keys():
            if oldindextonewindex[oendidx] not in atomswithcutbonds:
                atomswithcutbonds.append(oldindextonewindex[oendidx])
            continue
        if oendidx not in oldindextonewindex.keys() and obgnidx in oldindextonewindex.keys():
            if oldindextonewindex[obgnidx] not in atomswithcutbonds:
                atomswithcutbonds.append(oldindextonewindex[obgnidx])
            continue
        if oendidx not in oldindextonewindex.keys() and obgnidx not in oldindextonewindex.keys():
            continue
        endidx=oldindextonewindex[oendidx]
        bgnidx=oldindextonewindex[obgnidx]
        if isinring==True and endidx in aromaticindices and bgnidx in aromaticindices:
            aromaticbonds.append([endidx,bgnidx]) 
        bondorder=bond.GetBondType()
        mw.AddBond(bgnidx,endidx,bondorder)
    # STEP 8
    smarts=rdmolfiles.MolToSmarts(mw)
    # STEP 9
    for atom in mw.GetAtoms():
        atomidx=atom.GetIdx()
        if atomidx in aromaticindices:
            atom.SetIsAromatic(True)
    # STEP 10
    for bond in mw.GetBonds():
        endidx = bond.GetEndAtomIdx()
        bgnidx = bond.GetBeginAtomIdx()
        temp=[endidx,bgnidx]
        if temp in aromaticbonds or temp[::-1] in aromaticbonds:
            bond.SetIsAromatic(True)
    # STEP 11
    smartsfortransfer=rdmolfiles.MolToSmarts(mw)
    
    return smarts,smartsfortransfer


def GenerateFragmentSMARTSList(poltype,ls):
    """
    Intent: When wanting to find matches to amoeba09 database, need a way to generate many SMARTS of varying lengths centered around inpyt atom indices from input molecule. This is sort of like maximum common substructure search but only centered on atom desired so need custom code. 
    Input: List of atom indices from input molecule want to generate SMARTS strings for
    Output: List of SMARTS strings containing matches to part of molecule with input atom indices
    Referenced By: MatchAtomIndicesSMARTSToParameterSMARTS
    Description: 
    1. Generate all possible fragment indices containing input atom indices up to length 14 (4 SP3 atoms for torsion have max 14 atoms)
    2. Iterate over each list of atom indices just generated
    3. Generate a SMARTS string for each one and append to array 
    """
    fragsmartslist=[]
    # STEP 1
    atomindiceslist=GenerateAllPossibleFragmentIndices(poltype,ls,poltype.rdkitmol,14)
    # STEP 2
    for thels in atomindiceslist:
        # STEP 3
        smarts,smartsfortransfer=GenerateFragmentSMARTS(poltype,thels)
        if smartsfortransfer not in fragsmartslist:
            fragsmartslist.append(smartsfortransfer)


    return fragsmartslist


def FindHowManyBondTypes(poltype,trymol,atomidx):
    """
    Intent: Want to define a "type" that does not take into account global symmetry but only surrounding envioronment. Define it with atomicnumber, number of single,aromatic, double and triple bonds.
    Input: Mol object, atom index that want a type for
    Output: The atom type (array of atomic number, number of bond types for each bond type for that atom)
    Referenced By: MatchAllPossibleSMARTSToParameterSMARTS
    Description: 
    1. Grab the atom of interest from input atom index and mol object
    2. Iterate over neighbors of atom of interest
    3. Grab bond type and increase count of respective bond type (single, aromatic, double, triple)
    4. Construct atom type using information of atomic number as well as how many single, aromatic, double and triple bonds.
    """
    numsinglebonds=0
    numaromaticbonds=0
    numdoublebonds=0
    numtriplebonds=0
    # STEP 1
    atom=trymol.GetAtomWithIdx(atomidx)
    # STEP 2
    for natom in atom.GetNeighbors():
        natomidx=natom.GetIdx()
        bond=trymol.GetBondBetweenAtoms(atomidx,natomidx)
        # STEP 3
        bondtype=bond.GetBondTypeAsDouble() 
        if bondtype==1:
            numsinglebonds+=1
        elif bondtype==1.5:
            numaromaticbonds+=1
        elif bondtype==2:
            numdoublebonds+=1
        elif bondtype==3:
            numtriplebonds+=1
    atomicnum=atom.GetAtomicNum()
    # STEP 4
    atomtype=tuple([atomicnum,numsinglebonds,numaromaticbonds,numdoublebonds,numtriplebonds])

    return atomtype



def MatchAllPossibleSMARTSToParameterSMARTS(poltype,parametersmartslist,parametersmartstosmartslist,ls,mol,parametersmartstordkitmol,parametersmartstosmartsmatchingtoindices,parametersmartstomolsmatchingtoindices,parametersmartstordkitmolmatchingindices,parametersmartstoprmmolmatchingindices):
    """
    Intent: For a given input array of atom indices (one for atom, two for bond, three for angle etc) then try to match SMARTS containing those indices to amoeba09 database in order to find parameters suitable to the given environment.  
    Input: List of SMARTS strings from amoeba09, dictionary mapping amoeba09 SMARTS to final matched SMARTS from input molecule, input list of atom indices from molecule trying to find SMARTS matches to amoeba09 for, mol object, dictionary of amoeba09 SMARTS -> mol object, dictionary of amoeba09 SMARTS ->  matching indices of input molecule,  dictionary of amoeba09 SMARTS ->  matching indices of SMARTS generated from input molecule
    Output: dictionary mapping amoeba09 SMARTS to final matched SMARTS from input molecule, dictionary mapping amoeba09 SMARTS to boolean specifying if the SMARTS matches all neighbors of input list of atom indices in input molecule
    Referenced By: MatchAtomIndicesSMARTSToParameterSMARTS
    Description: 
    1. Iterate over amoeba09 SMARTS
    2. Grab amoeba09 SMARTS
    3. Grab mol object for amoeba09 SMARTS
    4. Grab list of SMARTS that matched to amoeba09 SMARTS
    5. Grab list of mols from SMARTS that matched to amoeba09 SMARTS
    6. Iterate over list of SMARTS that matched to amoeba09 SMARTS
    7. Grab SMARTS that matched to amoeba09 SMARTS
    8. Grab corresponding mol object for that SMARTS 
    9. Assume the SMARTS match is a "good" match for now
    10.Grab the first match of input molecule indices that matched from SMARTS 
    11.If the number of atoms in match is greater than or equal to number of atoms in input array of atom indices
    12.Iterate over all matches
    13.If there is an atom index in ls (input list of atom indices wish to find amoeba09 parameters for) that does not exist in match, then this is not a "good" match and skip that one. 
    14.Check if the matches are consecutively connected in the input molecule (otherwise SMARTS match wouldnt make sense)
    15.Define new array that includes neighbors of (ls) input array of atom indices. Determine if the match contains all the neighbors or not.
    16. Sanity check to ensure the number of rings in amoeba09 SMARTS doesnt exceed number of rings in SMARTS match from input mol (redundant now?)
    17. Grab the hybridizations of each atom in input molecule that matched to SMARTS, later will use to filter if match has same hybridzation as amoeba09 SMARTS
    18. Count the number of times each "type" defined by FindHowManyBondTypes, occurs for array of input atom indices + neighbors. Will later use as filtration step for potential SMARTS matches to amoeba09 parameter SMARTS
    19. Match generated SMARTS to input molecule
    20. Iterate over matches and if all indices in input array exist in match and match indices are consecutive
    21. Generate dictionary of input molecule index to SMARTS index (the atom matched in input molecule to index of atom in SMARTS)
    22. Grab indices that map current SMARTS to current amoeba09 SMARTS
    23. SMARTS matching from current SMARTS->amoeba09 SMARTS returns many matches, try and pick one that includes most neighbors from input array of atom indices from input molecule
    24. Skip any matches that dont have correct hybridization between input molecule atoms and atoms matched with amoeba09 SMARTS
    25. Compute various "scores" for filtering out possible SMARTS that match to amoeba09 SMARTS to find the "best" match. The first score is how many atoms from input array + neighbors are included in match (want to maximize this), the second score is the difference between number of local "types" in input molecule vs "types" in amoeba09 SMARTS (want to minimize this, a difference of 0 means molecule must be the same), the third score is the difference between number of each atomic number occuring in input molecule vs amoeba09 SMARTS (want to minimize).
    26. Apply filters sequentially.
    27. If no good matches are found in case of vdw matches, just try and match element only and keep as match (better than wild card match). 
    """
    smartsmcstomol={}
    prmsmartstomcssmarts={}
    parametersmartstoscore={}
    parametersmartstonumcommon={}
    parametersmartstootherscore={}
    parametersmartstothefinalscore={}
    parametersmartstofoundallneighbs={}
    newls=copy.deepcopy(ls)
    # STEP 1
    for parametersmartsidx in tqdm(range(len(parametersmartslist)),desc='amoeba09 search for '+str(ls)):
        # STEP 2
        parametersmarts=parametersmartslist[parametersmartsidx]
        # STEP 3
        prmmol=parametersmartstordkitmol[parametersmarts]
        # STEP 4
        smartsmatchingtoindices=parametersmartstosmartsmatchingtoindices[parametersmarts]
        # STEP 5
        molsmatchingtoindices=parametersmartstomolsmatchingtoindices[parametersmarts]
        prmsmartsatomnum=prmmol.GetNumAtoms()
        thesmartstonumneighbs={}
        thesmartstotypescore={}
        thesmartstoelementscore={}
        thesmartstothemol={}
        # STEP 6
        for theidx in range(len(smartsmatchingtoindices)):
            # STEP 7
            thesmarts=smartsmatchingtoindices[theidx]
            # STEP 8
            themol=molsmatchingtoindices[theidx]
            thesmartstothemol[thesmarts]=themol
            # STEP 9
            diditmatchprmmol=True
            diditmatch=True
            if diditmatch==True and diditmatchprmmol==True:
                matches=parametersmartstordkitmolmatchingindices[parametersmarts][theidx]
                # STEP 10
                firstmatch=matches[0]
                # STEP 11 
                if len(firstmatch)>=len(ls):
                    # STEP 12
                    for match in matches:
                        goodmatch=True
                        matchidxs=[]
                        # STEP 13
                        for idx in ls:
                            if idx not in match:
                                goodmatch=False
                                break
                            else:
                                matchidx=match.index(idx)
                                matchidxs.append(matchidx)
                        # STEP 14
                        if len(ls)>1 and goodmatch==True:
                            goodmatch=CheckIfConsecutivelyConnected(poltype,matchidxs,themol)
                        # STEP 15
                        newls=ExtendByNeighbors(poltype,ls) 
                        score=0 
                        allneighbsin=True
                        matchidxs=[]
                        for idx in newls:
                            if idx in match:
                                score+=1
                                matchidx=match.index(idx)
                                matchidxs.append(matchidx)
                            else:
                                allneighbsin=False
                        if goodmatch==True:
                             
                            break
                    # STEP 16 
                    prmmolnumrings=CountRingsInSMARTS(poltype,parametersmarts)
                    molnumrings=CalcNumRings(poltype.rdkitmol)
                    if prmmolnumrings>molnumrings:
                        goodmatch=False


                    if goodmatch==True:
                        # STEP 17
                        rdkitatoms=[poltype.rdkitmol.GetAtomWithIdx(r) for r in firstmatch]
                        rdkithybs=[r.GetHybridization() for r in rdkitatoms]
                        rdkitmatch=copy.deepcopy(firstmatch)
                        idx=ls[0]
                            


                        # STEP 18
                        rdkitatomictypetotype={}
                        rdkitatomicnumtonum={}
                        for atom in poltype.rdkitmol.GetAtoms():
                            atomicnum=atom.GetAtomicNum()
                            atomidx=atom.GetIdx()
                            if atomicnum not in rdkitatomicnumtonum.keys():
                                rdkitatomicnumtonum[atomicnum]=0
                            rdkitatomicnumtonum[atomicnum]+=1
                            if atomidx in newls:
                                atomtype=FindHowManyBondTypes(poltype,poltype.rdkitmol,atomidx)
                                
                                if atomtype not in rdkitatomictypetotype.keys():
                                    rdkitatomictypetotype[atomtype]=0
                                rdkitatomictypetotype[atomtype]+=1
                        # STEP 19 
                        sp=openbabel.OBSmartsPattern()
                        openbabel.OBSmartsPattern.Init(sp,thesmarts)
                        diditmatch=sp.Match(poltype.mol)
                        babelmatches=sp.GetMapList()
                        newbabelmatches=[]
                        for mtch in babelmatches:
                            newmtch=[i-1 for i in mtch]
                            newbabelmatches.append(newmtch)
                        for newmtch in newbabelmatches:
                            if newmtch not in matches:
                                matches.append(newmtch)
                        # STEP 20
                        for match in matches:
                            validmatch=CheckMatch(poltype,match,ls,thesmarts,themol)
                            if validmatch==True:
                               indices=list(range(len(match)))
                               smartsindextomoleculeindex=dict(zip(indices,match)) 
                               moleculeindextosmartsindex={v: k for k, v in smartsindextomoleculeindex.items()}


                        # STEP 21
                        smartindices=[moleculeindextosmartsindex[i] for i in ls]
                        # STEP 22
                        prmmatches=parametersmartstoprmmolmatchingindices[parametersmarts][theidx]
                        # STEP 23 
                        firstmatch=TryAndPickMatchWithNeighbors(poltype,prmmatches,smartindices,themol)
                        # STEP 24
                        indices=list(range(len(firstmatch)))
                        smartsindextoparametersmartsindex=dict(zip(indices,firstmatch)) 
                        prmsmartsindices=[smartsindextoparametersmartsindex[i] for i in smartindices]
                        prmatoms=[prmmol.GetAtomWithIdx(r) for r in firstmatch]
                        prmhybs=[r.GetHybridization() for r in prmatoms]
                        allhybssame=True
                        for r in range(len(rdkithybs)):
                            rdkithyb=rdkithybs[r]
                            prmhyb=prmhybs[r]
                            if rdkithyb!=prmhyb:
                                allhybssame=False
                        if allhybssame==False:
                            continue
                        # STEP 25 
                        prmsmartsatomicnumtonum={}
                        for atomicnum,num in rdkitatomicnumtonum.items():
                            if atomicnum not in prmsmartsatomicnumtonum.keys():
                                prmsmartsatomicnumtonum[atomicnum]=0
                        prmsmartsatomictypetotype={}
                        for atomicnum,num in rdkitatomictypetotype.items():
                            prmsmartsatomictypetotype[atomicnum]=0
                        extra=prmsmartsindices[:]
                        for atom in prmmol.GetAtoms():
                            atomicnum=atom.GetAtomicNum()
                            atomidx=atom.GetIdx()
                            if atomicnum not in prmsmartsatomicnumtonum.keys():
                                prmsmartsatomicnumtonum[atomicnum]=0
                            prmsmartsatomicnumtonum[atomicnum]+=1
                            if atomidx in prmsmartsindices:
                                atomtype=FindHowManyBondTypes(poltype,prmmol,atomidx)

                                if atomtype not in prmsmartsatomictypetotype.keys():
                                    prmsmartsatomictypetotype[atomtype]=0
                                prmsmartsatomictypetotype[atomtype]+=1
                                for natom in atom.GetNeighbors():
                                    natomidx=natom.GetIdx()
                                    atomicnum=natom.GetAtomicNum()

                                    if natomidx not in extra:
                                        atomtype=FindHowManyBondTypes(poltype,prmmol,natomidx)
                                        if atomtype not in prmsmartsatomictypetotype.keys():
                                            prmsmartsatomictypetotype[atomtype]=0
                                        prmsmartsatomictypetotype[atomtype]+=1
                                        extra.append(natomidx)
                                    if len(ls)==1 and len(atom.GetNeighbors())==1:
                                        for nnatom in natom.GetNeighbors():
                                            nnatomidx=nnatom.GetIdx()
                                            atomicnum=nnatom.GetAtomicNum()

                                            if nnatomidx not in extra:
                                                atomtype=FindHowManyBondTypes(poltype,prmmol,nnatomidx)
                                                if atomtype not in prmsmartsatomictypetotype.keys():
                                                    prmsmartsatomictypetotype[atomtype]=0
                                                prmsmartsatomictypetotype[atomtype]+=1
                                                extra.append(natomidx)




                        otherscore=0
                        for atomicnum,num in rdkitatomictypetotype.items():
                            prmnum=prmsmartsatomictypetotype[atomicnum]
                            diff=np.abs(prmnum-num)

                            otherscore+=diff

                        globalscore=0
                        for atomicnum,num in rdkitatomicnumtonum.items():
                            prmnum=prmsmartsatomicnumtonum[atomicnum]
                            diff=np.abs(prmnum-num)

                            globalscore+=diff

                        thesmartstonumneighbs[thesmarts]=score
                        thesmartstotypescore[thesmarts]=otherscore
                        thesmartstoelementscore[thesmarts]=globalscore
    

        if len(thesmartstonumneighbs.keys())>0:
            maxscore=max(thesmartstonumneighbs.values())            
            for thesmarts,score in thesmartstonumneighbs.items():
                themol=thesmartstothemol[thesmarts]
                otherscore=thesmartstotypescore[thesmarts]
                if score==maxscore:

                    prmsmartstomcssmarts[parametersmarts]=thesmarts
                    parametersmartstoscore[parametersmarts]=score
                    parametersmartstootherscore[parametersmarts]=otherscore
                    smartsmcstomol[thesmarts]=themol
                    parametersmartstofoundallneighbs[parametersmarts]=allneighbsin
                    parametersmartstothefinalscore[parametersmarts]=thesmartstoelementscore[thesmarts]
                    break
    foundmin=False
    # STEP 26
    parametersmartstofinalscore={}
    parametersmartstolastscore={}
    if len(parametersmartstoscore.keys())>0:
        minscore=max(parametersmartstoscore.values())
        for parametersmarts in parametersmartstoscore.keys():
            score=parametersmartstoscore[parametersmarts]
            otherscore=parametersmartstootherscore[parametersmarts]
            smartsmcs=prmsmartstomcssmarts[parametersmarts]
            lastscore=parametersmartstothefinalscore[parametersmarts]
            if score==minscore:
                parametersmartstolastscore[parametersmarts]=lastscore          
                parametersmartstofinalscore[parametersmarts]=otherscore
        minscore=min(parametersmartstofinalscore.values())
        finalprmsmartstoscore={}
        for parametersmarts in parametersmartstofinalscore.keys():
            score=parametersmartstofinalscore[parametersmarts]
            smartsmcs=prmsmartstomcssmarts[parametersmarts]
            lastscore=parametersmartstolastscore[parametersmarts]
            if score==minscore:
                finalprmsmartstoscore[parametersmarts]=lastscore
                
        minscore=min(finalprmsmartstoscore.values())
        for parametersmarts,score in finalprmsmartstoscore.items():
            if score==minscore:
                smartsmcs=prmsmartstomcssmarts[parametersmarts] 
                mcsmol=smartsmcstomol[smartsmcs]
                foundmin=True
                break

        if foundmin==True:
            smartls=[smartsmcs,smartsmcs]
            if parametersmarts not in parametersmartstosmartslist.keys():
                parametersmartstosmartslist[parametersmarts]=smartls
    else:
        # STEP 27
        if len(ls)==1:
            atom=poltype.rdkitmol.GetAtomWithIdx(ls[0])
            atomicnum=atom.GetAtomicNum()
            string='[#'+str(atomicnum)+']'
            othermol=Chem.MolFromSmarts(string)
            for parametersmarts in parametersmartslist:
                prmmol=Chem.MolFromSmarts(parametersmarts)
                mols = [othermol,prmmol]
                diditmatch=mols[1].HasSubstructMatch(mols[0])
                if diditmatch==True:
                    matches=mols[1].GetSubstructMatches(mols[0])
                    firstmatch=matches[0]
                    smartls=[string,string]
                    if parametersmarts not in parametersmartstosmartslist.keys():
                        parametersmartstosmartslist[parametersmarts]=smartls
                        parametersmartstofoundallneighbs[parametersmarts]=False
                    break

         

    return parametersmartstosmartslist,parametersmartstofoundallneighbs


def GenerateAllPossibleFragmentIndices(poltype,ls,rdkitmol,maxatomsize):
    """
    Intent: When generating SMARTS from input molecule to match to amoeba09 SMARTS, need to generate many atomic indices first for each possible SMARTS.
    Input: Input array of atom indices from input molecule, mol object, maximum number of atoms to search out in molecule.
    Output: List of list of atomic indices each containing input array of atom indices
    Referenced By: GenerateFragmentSMARTSList
    Description: 
    1. Start with list of current atom indices
    2. Define a second list that will constantly be updated with more atoms and stop when the lists are the same. After first iteration set oldindexlist=indexlist, if at end of iteration, no new neighbors are added, then will exit while loop.
    3. If length of list is greater than maxatomsize or if searched further than two neighbors away, then quit
    4. Grab neighboring indices and for each neighbor, generate a list of original atom indices+some of neigbors, append to list of atom indices to generate SMARTS from
    5. Add neighboring indices to indexlist 
    """
    # STEP 1
    atomindiceslist=[copy.deepcopy(ls)]
    oldindexlist=copy.deepcopy(ls)
    indexlist=[]
    count=0
    neighbcount=0
    # STEP 2
    while set(oldindexlist)!=set(indexlist):
        if count!=0:
            oldindexlist=copy.deepcopy(indexlist)
        # STEP 3
        if len(oldindexlist)>maxatomsize:
            break
        if neighbcount>=2:
            break
        # STEP 4
        neighborindexes=GrabNeighboringIndexes(poltype,oldindexlist,rdkitmol)
        neighbcount+=1
       
        for i in range(len(neighborindexes)):
            combs=list(itertools.combinations(neighborindexes, i+1))
            for comb in combs:
                newcomb=copy.deepcopy(oldindexlist)
                newcomb.extend(comb) 
                if newcomb not in atomindiceslist:
                    atomindiceslist.append(newcomb)
        count+=1
        # STEP 5
        newindexlist=copy.deepcopy(oldindexlist)
        for index in neighborindexes:
            if index not in newindexlist:
                newindexlist.append(index)
        indexlist=newindexlist
  
    return atomindiceslist


def GrabNeighboringIndexes(poltype,indexlist,rdkitmol):
    """
    Intent: Grab neighboring indices for when generating all possible SMARTS to match to a trial amoeba09 SMARTS.
    Input: Array of atom indices, mol object
    Output: Array of neighboring atom indices
    Referenced By: GenerateAllPossibleFragmentIndices
    Description: 
    1. Iterate over atom indices in input list
    2. Grab neighbors of each atom
    3. Append to array of neighbors 
    """
    neighborindexes=[]
    # STEP 1
    for index in indexlist:
        atom=rdkitmol.GetAtomWithIdx(index)
        # STEP 2
        for natom in atom.GetNeighbors():
            natomidx=natom.GetIdx()
            # STEP 3
            if natomidx not in neighborindexes and natomidx not in indexlist:
                neighborindexes.append(natomidx)
    return neighborindexes


def CountRingsInSMARTS(poltype,parametersmarts):
    """
    Intent: Sanity check to ensure number of rings in amoeba09 SMARTS, doesnt exceed number of rings in input molecule 
    Input: SMARTS strings
    Output: Number of rings
    Referenced By: MatchAllPossibleSMARTSToParameterSMARTS
    Description: 
    1. Iterate over SMARTS characters
    2. Search for when after ] bracket (after atom) and if there is a digit after that (specifies ring number), append to array of ring numbers
    3. Return length of ring numbers array 
    """
    closedbrack=True
    numbers=[]
    # STEP 1
    for e in parametersmarts:
        # STEP 2
        if e=='[':
            closedbrack=False
        elif e==']':
            closedbrack=True
        elif e.isdigit():
            if closedbrack==True:
                if e not in numbers:
                    numbers.append(e)
    # STEP 3
    rings=len(numbers)
    return rings

def CheckIfConsecutivelyConnected(poltype,matchidxs,mcsmol):
    """
    Intent: Check if indices in match are consecutively connected (needed for adequate SMARTS match to amoeba09 SMARTS).
    Input: Array of indices in SMARTS that match to input molecule, mol from SMARTS 
    Output: Boolean specifying if all indices are consecutive or not.
    Referenced By: MatchAllPossibleSMARTSToParameterSMARTS
    Description: 
    1. Iterate over atom indices (except last one)
    2. Grab atom from SMARTS mol that corresponds to index
    3. Grab next index
    4. Grab neighboring indices of current index
    5. If the next atom index in sequence does not exist in neighboring indices of current index, then the match is not consecutively connected. 
    """
    goodmatch=True
    # STEP 1
    for i in range(len(matchidxs)-1):
        matchidx=matchidxs[i]
        # STEP 2
        atom=mcsmol.GetAtomWithIdx(matchidx)
        # STEP 3
        nextmatchidx=matchidxs[i+1]
        natoms=atom.GetNeighbors()
        # STEP 4
        natomidxs=[a.GetIdx() for a in natoms]
        # STEP 5
        if nextmatchidx not in natomidxs:
            goodmatch=False


    return goodmatch





def MatchAtomIndicesSMARTSToParameterSMARTS(poltype,listforprm,parametersmartslist,mol,parametersmartstordkitmol):
    """
    Intent: For each list of atom indices (all atoms, all bonds, all angles etc), find the "best" amoeba09 SMARTS match. 
    Input: List of list of atom indices, list of all amoeba09 SMARTS, mol object, dictionary of amoeba09 SMARTS to corresponding mol object.
    Output: Dictionary of atom indices (vdw, bond, angle, torsion) -> best amoeba09 SMARTS match, Dictionary of atom indices (vdw, bond, angle, torsion) -> boolean if amoeba09 SMARTS matches to all neighbors of atom indices 
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over each list of atom indices
    2. For each list of atom indices, generate a list of potential SMARTS to match to any amoeba09 SMARTS
    3. Save results in list for later 
    4. Iterate over list of amoeba09 SMARTS
    5. Iterate over list of all possible SMARTS to match to current amoeba09 SMARTS
    6. If the current SMARTS matches to both input molecule and to amoeba09 SMARTS molecule, then will save for later filtering
    7. Save SMARTS match results in dictionaries for later processing and filtering
    8. Iterate over each list of atom indices again
    9. Call MatchAllPossibleSMARTSToParameterSMARTS using current list of atom indices and the saved list of amoeba09 SMARTS and potential SMARTS matches from input molecule. This will output the "best" SMARTS match for each amoeba09 SMARTS.
    10. Choose the amoeba09 SMARTS, that has the largest SMARTS matching to it (most information in input molecule environment)
    11. If no match was found, assign a wildcard SMARTS (worst case)  
    """
    listforprmtoparametersmarts={}
    listforprmtosmarts={}
    listforprmtomatchallneighbs={}
    allfragsmartslist=[]
    # STEP 1
    for ls in listforprm:
        # STEP 2
        fragsmartslist=GenerateFragmentSMARTSList(poltype,ls)
        # STEP 3
        for fragsmarts in fragsmartslist:
            if fragsmarts not in allfragsmartslist:
                allfragsmartslist.append(fragsmarts)

    parametersmartstosmartsmatchingtoindices={}
    parametersmartstomolsmatchingtoindices={}
    parametersmartstordkitmolmatchingindices={}
    parametersmartstoprmmolmatchingindices={}
    # STEP 4
    for parametersmarts in parametersmartslist:
        prmmol=parametersmartstordkitmol[parametersmarts]
        finalfragsmartslist=[]
        smartsmatchingtoindices=[]
        molsmatchingtoindices=[]
        rdkitmolmatchingindices=[]
        prmmolmatchingindices=[]
        # STEP 5
        for fragsmarts in allfragsmartslist:
            fragmol=Chem.MolFromSmarts(fragsmarts)
            diditmatch=poltype.rdkitmol.HasSubstructMatch(fragmol)
            # STEP 6
            if diditmatch==True:
                diditmatch=prmmol.HasSubstructMatch(fragmol)
                if diditmatch==True:
                    matches=list(poltype.rdkitmol.GetSubstructMatches(fragmol,maxMatches=10000))
                    prmmatches=prmmol.GetSubstructMatches(fragmol)
                    molsmatchingtoindices.append(fragmol)
                    smartsmatchingtoindices.append(fragsmarts)
                    rdkitmolmatchingindices.append(matches)
                    prmmolmatchingindices.append(prmmatches)
        # STEP 7
        parametersmartstosmartsmatchingtoindices[parametersmarts]=smartsmatchingtoindices
        parametersmartstomolsmatchingtoindices[parametersmarts]=molsmatchingtoindices
        parametersmartstordkitmolmatchingindices[parametersmarts]=rdkitmolmatchingindices
        parametersmartstoprmmolmatchingindices[parametersmarts]=prmmolmatchingindices
    # STEP 8
    for ls in listforprm:
        # STEP 9 
        parametersmartstomatchlen={}
        parametersmartstosmartslist={}
        parametersmartstosmartslist,parametersmartstofoundallneighbs=MatchAllPossibleSMARTSToParameterSMARTS(poltype,parametersmartslist,parametersmartstosmartslist,ls,mol,parametersmartstordkitmol,parametersmartstosmartsmatchingtoindices,parametersmartstomolsmatchingtoindices,parametersmartstordkitmolmatchingindices,parametersmartstoprmmolmatchingindices)
        # STEP 10
        if len(parametersmartstosmartslist.keys())!=0:
            parametersmartstosmartslen={}
            for prmsmarts,newls in parametersmartstosmartslist.items():
                smarts=newls[0]
                smartslen=len(smarts)
                parametersmartstosmartslen[prmsmarts]=smartslen
            
            valuelist=list(parametersmartstosmartslen.values())
            maxvalue=max(valuelist)
            for prmsmarts,smartslen in parametersmartstosmartslen.items():
                if smartslen==maxvalue:
                    maxprmsmarts=prmsmarts
                    break 
            maxsmartsls=parametersmartstosmartslist[maxprmsmarts]
            matchallneighbs=parametersmartstofoundallneighbs[maxprmsmarts]
        else:
            # STEP 11
            matchallneighbs=False
            wild='[*]'
            smarts=''
            for i in range(len(ls)):
                smarts+=wild+'~'
            smarts=smarts[:-1]
            maxsmartsls=[smarts,smarts]
            maxprmsmarts='[#6](-[H])(-[H])(-[H])-[#6](-[H])(-[H])(-[H])'
        
        listforprmtoparametersmarts[tuple(ls)]=maxprmsmarts
        listforprmtosmarts[tuple(ls)]=maxsmartsls
        listforprmtomatchallneighbs[tuple(ls)]=matchallneighbs
    return listforprmtoparametersmarts,listforprmtosmarts,listforprmtomatchallneighbs



def CheckMatch(poltype,match,atomindices,smarts,substructure):
    """
    Intent: Want to ensure that when matching SMARTS to input molecule, the atom indices list (vdw, bond, angle etc) are all consecutively matched in the SMARTS match. 
    Input: Array of indices in input molecule that correspond to the SMARTS match, array of indices of interest from input molecule, SMARTS string, mol for SMARTS string 
    Output: Boolean specifying if match is consecutive or not
    Referenced By: MatchAllPossibleSMARTSToParameterSMARTS , GenerateAtomIndexToAtomTypeAndClassForAtomList
    Description: 
    1. If all indices in molecule of interest are not in match array, then its not a valid match to begin with
    2. Assume indices are consecutive
    3. For each index in indices of interest, find location of that index in match array
    4. Save each location found in new array
    5. Now using the mol and location indices array, check if each atom that should be consecutive is consecutive in that match array, this will return boolean if consecutive or not
    """
    # STEP 1
    allin=True
    for idx in atomindices:
        if idx not in match:
            allin=False
    # STEP 2
    validmatch=True
    if allin==True:
        indices=[]
        # STEP 3
        for idx in atomindices:
            matchidx=match.index(idx)
            # STEP 4
            indices.append(matchidx)
        # STEP 5
        checkconsec=CheckConsecutiveTorsion(poltype,indices,substructure)
        if checkconsec==False:
            validmatch=False
    else:
        validmatch=False
    return validmatch 
       
def CheckConsecutiveTorsion(poltype,indices,substructure):
    """
    Intent: Check for each parameter type, bond, angle, torsion, if the match to SMARTS is consecutive. 
    Input: Array of atom indices of interest in molecule sorted by location found in SMARTS match, mol for SMARTS
    Output: Boolean if consecutive or not
    Referenced By: CheckMatch
    Description: 
    1. For bond, just check if there exists bond between indices or not
    2. For angle, just check if there exists bonds between a,b and b,c or not
    3. For torsion, just check if there exists bonds between a,b and b,c and c,d or not
    """
    consec=True
    # STEP 1
    if len(indices)==2:
        a,b=indices[:]
        bond=substructure.GetBondBetweenAtoms(a,b)
        if bond==None:
            consec=False
    # STEP 2
    elif len(indices)==3:
        a,b,c=indices[:]
        bond=substructure.GetBondBetweenAtoms(a,b)
        if bond==None:
            consec=False
        bond=substructure.GetBondBetweenAtoms(b,c)
        if bond==None:
            consec=False
    # STEP 3
    elif len(indices)==4: 
        a,b,c,d=indices[:]
        bond=substructure.GetBondBetweenAtoms(a,b)
        if bond==None:
            consec=False
        bond=substructure.GetBondBetweenAtoms(b,c)
        if bond==None:
            consec=False
        bond=substructure.GetBondBetweenAtoms(c,d)
        if bond==None:
            consec=False


    return consec


def GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol):
    """
    Intent: For each list of atom indices with the "best" amoeba09 SMARTS match, determine the tinker type and class number and store in dictionaries for later use.
    Input: Dictionary of atom indices -> "best" amoeba09 SMARTS ,  Dictionary of atom indices -> "best" SMARTS that matches to amoeba09 SMARTS and input molecule, dictionary of amoeba09 SMARTS, element -> tinker type description, dictionary of element tinker description -> tinker type, dictionary of tinker type -> tinker class, input mol object
    Output: Dictionary of atom indices -> tinker types, Dictionary of atom indices -> tinker classes, Dictionary of atom indices -> amoeba09 SMARTS + atom order, Dictionary of atom indices -> tinker description + element, Dictionary of atom indices -> SMARTS + atom order
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over dictionary of list of atom indices to best amoeba09 SMARTS
    2. Grab corresponding best SMARTS match to input molecule 
    3. Generate mol object from SMARTS
    4. Generate matches between SMARTS and input molecule
    5. Iterate over matches
    6. Check if match is valid (consecutively matches atom indices in molecule).
    7. If match is valid, then generate a dictionary of input molecule index -> index of match in SMARTS 
    8. Generate rdkit mol object from amoeba09 SMARTS
    9. Output .mol file and have babel read as input to generate babel mol object
    10. Compute dictionary of atom index (in amoeba09 SMARTS) -> type number, will later use to help matching indices of same atom type 
    11. Generate matches for SMARTS to amoeba09 SMARTS
    12. Try and pick match from 11. that includes neighbors of atom indices of interest in SMARTS
    13. Generate dictionary of SMARTS atom index -> amoeba09 SMARTS atom index
    14. Generate a dictionary of amoeba09 SMARTS and atom order (location of atoms of interest in amoeba09 SMARTS) -> element + tinker description using input dictionaries and current amoeba09 SMARTS
    15. Using the type dictionary generated earlier, iterate over amoeba09 molecule and for atoms of same type, fill in dictionary from 14 with information of atoms having same type
    16. Using dictionary from 15 and input dictionaries, fill in desired dictionaries mapping atom indices -> tinker types/classes/amoeba09 SMARTS/SMARTS  
    """
    atomindicestotinkertypes={}
    atomindicestotinkerclasses={}
    atomindicestoparametersmartsatomorders={}
    atomindicestoelementtinkerdescrips={}
    atomindicestosmartsatomorders={}
    # STEP 1
    for atomindices,parametersmarts in atomindicesforprmtoparametersmarts.items():
        # STEP 2
        smartsls=atomindicesforprmtosmarts[atomindices]
        smarts=smartsls[0]
        smartsfortransfer=smartsls[1]
        # STEP 3
        substructure = Chem.MolFromSmarts(smarts)
        # STEP 4
        matches=list(rdkitmol.GetSubstructMatches(substructure,maxMatches=10000))
        sp=openbabel.OBSmartsPattern()
        openbabel.OBSmartsPattern.Init(sp,smarts)
        diditmatch=sp.Match(poltype.mol)
        babelmatches=sp.GetMapList()
        newbabelmatches=[]
        for mtch in babelmatches:
            newmtch=[i-1 for i in mtch]
            newbabelmatches.append(newmtch)
        for newmtch in newbabelmatches:
            if newmtch not in matches:
                matches.append(newmtch)
        # STEP 5
        for match in matches:
            # STEP 6
            validmatch=CheckMatch(poltype,match,atomindices,smarts,substructure)
            # STEP 7
            if validmatch==True:
               indices=list(range(len(match)))
               smartsindextomoleculeindex=dict(zip(indices,match)) 
               moleculeindextosmartsindex={v: k for k, v in smartsindextomoleculeindex.items()}
        # STEP 8
        structure = Chem.MolFromSmarts(parametersmarts)
        fragmentfilepath='fragment.mol'
        if os.path.isfile(fragmentfilepath):
            os.remove(fragmentfilepath)
        # STEP 9
        rdmolfiles.MolToMolFile(structure,fragmentfilepath)
        obConversion = openbabel.OBConversion()
        fragbabelmol = openbabel.OBMol()
        inFormat = obConversion.FormatFromExt(fragmentfilepath)
        obConversion.SetInFormat(inFormat)
        obConversion.ReadFile(fragbabelmol, fragmentfilepath)
        # STEP 10
        fragidxtosymclass,symmetryclass=symm.gen_canonicallabels(poltype,fragbabelmol,rdkitmol=structure)
        smartindices=[moleculeindextosmartsindex[i] for i in atomindices]
        # STEP 11
        substructure = Chem.MolFromSmarts(smartsfortransfer)
        matches=structure.GetSubstructMatches(substructure)
        # STEP 12
        firstmatch=TryAndPickMatchWithNeighbors(poltype,matches,smartindices,substructure)
        indices=list(range(len(firstmatch)))
        # STEP 13
        smartsindextoparametersmartsindex=dict(zip(indices,firstmatch)) 
        # STEP 14
        parametersmartsordertoelementtinkerdescrip={}
        for parametersmartsatomorderlist,elementtinkerdescrip in smartsatomordertoelementtinkerdescrip.items():
            prmsmarts=parametersmartsatomorderlist[0]
            atomorderlist=parametersmartsatomorderlist[1]
            if prmsmarts==parametersmarts:
                atomorderlist=parametersmartsatomorderlist[1]
                for atomorder in atomorderlist:
                    parametersmartsordertoelementtinkerdescrip[atomorder]=elementtinkerdescrip 
        # STEP 15
        for fragidx,symclass in fragidxtosymclass.items():
            indexes=GrabKeysFromValue(poltype,fragidxtosymclass,symclass)
            specialindex=None
            for index in indexes:
                if index in parametersmartsordertoelementtinkerdescrip.keys():
                    specialindex=index
            if specialindex!=None:
                elementtinkerdescrip=parametersmartsordertoelementtinkerdescrip[specialindex]
                for index in indexes:
                    parametersmartsordertoelementtinkerdescrip[index]=elementtinkerdescrip
        # STEP 16
        parametersmartindices=[smartsindextoparametersmartsindex[i] for i in smartindices]
        parametersmartsorders=[i+1 for i in parametersmartindices]
        elementtinkerdescrips=[parametersmartsordertoelementtinkerdescrip[i] for i in parametersmartsorders]
        tinkertypes=[elementtinkerdescriptotinkertype[i] for i in elementtinkerdescrips]
        tinkerclasses=[tinkertypetoclass[i] for i in tinkertypes]
        tinkerclasses=[int(i) for i in tinkerclasses]
        tinkertypes=[int(i) for i in tinkertypes]
        smartsorders=[i+1 for i in smartindices]
        atomindicestotinkertypes[atomindices]=tinkertypes
        atomindicestotinkerclasses[atomindices]=tinkerclasses
        atomindicestoparametersmartsatomorders[atomindices]=[parametersmarts,parametersmartsorders]
        atomindicestoelementtinkerdescrips[atomindices]=elementtinkerdescrips             
        atomindicestosmartsatomorders[atomindices]=[smarts,smartsorders]

    return atomindicestotinkertypes,atomindicestotinkerclasses,atomindicestoparametersmartsatomorders,atomindicestoelementtinkerdescrips,atomindicestosmartsatomorders

def TryAndPickMatchWithNeighbors(poltype,matches,smartindices,substructure):
    """
    Intent: Trying to pick a match from SMARTS -> amoeba09 SMARTS that includes more neighbors (corresponds to including more neighbors of input molecule around atom indices of interest)
    Input: List of list of atom indices that match to amoeba09 SMARTS from SMARTS, the indices in SMARTS that correspond to same indices of interst in input molecule 
    Output: Match that attempt to include neighbors of SMARTS mol matches to indices of interest in input molecule 
    Referenced By: MatchAllPossibleSMARTSToParameterSMARTS,GenerateAtomIndexToAtomTypeAndClassForAtomList
    Description:
    1. Generate array of atom indices from SMARTS mol that include neighbors of input smartindices
    2. Iterate over each match and check if all indices from step 1. are included in match. If so, then use this match. 
    """
    # STEP 1
    extendedsmartindices=[]
    for idx in smartindices:
        atom=substructure.GetAtomWithIdx(idx)
        extendedsmartindices.append(idx)
        for natm in atom.GetNeighbors():
            natmidx=natm.GetIdx()
            if natmidx not in extendedsmartindices:
                extendedsmartindices.append(natmidx)
    # STEP 2
    for match in matches:
        indices=list(range(len(match)))
        smartsindextoparametersmartsindex=dict(zip(indices,match)) 
        foundmap=True
        for idx in extendedsmartindices:
            if idx not in smartsindextoparametersmartsindex:
                foundmap=False
        if foundmap==True:
            return match
    firstmatch=matches[0]
    return firstmatch
   
        

def GrabSMARTSList(poltype,smartsatomordertoelementtinkerdescrip):
    """
    Intent: Just extract SMARTS list from dictionary containing SMARTS + atom order -> element + tinker description
    Input: dictionary containing SMARTS + atom order -> element + tinker description
    Output: list of SMARTS
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over input dictionary
    2. Extract SMARTS from key
    3. Save to list 
    """
    smartslist=[]
    # STEP 1
    for smartsatomorder in smartsatomordertoelementtinkerdescrip.keys():
        # STEP 2
        smarts=smartsatomorder[0]
        # STEP 3
        if smarts not in smartslist:
            smartslist.append(smarts)
    return smartslist 

def CheckForPlanerAngles(poltype,listofanglesforprm,mol):
    """
    Intent: Tinker requires angle parameters that should be planer to be anglep rather than just angle, so need to detect this.
    Input: List of all angles in input molecule, mol object
    Output: List of angles in input molecule that are planar
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over each angle in molecule
    2. Check if the hybridization of the middle atom is SP2, if so and that atom has 3 neighbors, then this is a planar angle. 
    """
    listofanglesthatneedplanarkeyword=[]
    # STEP 1
    for ls in listofanglesforprm:
        a = mol.GetAtom(ls[0]+1)
        b = mol.GetAtom(ls[1]+1)
        c = mol.GetAtom(ls[2]+1)
        # STEP 2
        anglep=False
        if b.GetHyb()==2: # only for SP2 hyb middle atoms use angp
            neighbs=list(openbabel.OBAtomAtomIter(b))
            if len(neighbs)==3:
                anglep=True
        if anglep==True:
            listofanglesthatneedplanarkeyword.append(ls)
    return listofanglesthatneedplanarkeyword

def ModifyAngleKeywords(poltype,angleprms,listofanglesthatneedplanarkeywordtinkerclassestopoltypeclasses):
    """
    Intent: Modify angle parameters taken from prm file to have anglep keyword if the angle is planer in input molecule
    Input: Array of angle parameter lines, dictionary of planar angles tinker classes -> poltype classes 
    Output: Modified array of angle parameter lines
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over each angle parameter line
    2. Grab array of tinker classes from parameter line
    3. Iterate over input dictionary
    4. If there is a match between tinker classes in parameter line and current iteration over dictionary, then replace angle with anglep 
    5. Append modified parameter line to array
    """
    newangleprms=[]
    # STEP 1
    for line in angleprms:
        found=False
        linesplit=line.split()
        # STEP 2
        temp=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
        # STEP 3
        for ls,polclassesls in listofanglesthatneedplanarkeywordtinkerclassestopoltypeclasses.items():
            # STEP 4
            inline=True
            for i in temp:
                if i not in polclassesls:
                    inline=False
            if inline==True:
                found=True
                if 'anglep' not in line:
                    newline=line.replace('angle','anglep')
                else:
                    newline=line
                break
        newline=line.replace('anglef','angle')

        if found==False:
            newline=newline.replace('anglep','angle')
        # STEP 5
        newangleprms.append(newline)
    return newangleprms

def FilterList(poltype,allitems,listbabel):
    """
    Intent: When wanting to replace bond length/angle lengths with QM optimized bond length/angle lengths, need a way to grab all bonds of same type but also ensure they are physically connected ( you can enumerate many combinations of indices for each bond type but not all are physically connected). 
    Input: All possible bonds/angles, list of physical bonds/angles
    Output: List of bonds/angles with same types that are also in input molecule
    Referenced By: AddOptimizedBondLengths , AddOptimizedAngleLengths
    Description: 
    1. Iterate over each bond/angle in input list
    2. Check if the bond/angle exists in molecule, if so append to array 
    """
    newallitems=[]
    # STEP 1
    for ls in allitems:
        revls=ls[::-1]
        # STEP 2
        if list(ls) in listbabel or list(revls) in listbabel: # need to check reverse too
            newallitems.append(ls)
    return newallitems

def AddOptimizedBondLengths(poltype,optmol,bondprms,bondlistbabel):
    """
    Intent: Modify bond parameters to use QM geometry optimized bond length values 
    Input: QM optimized mol object, list of bond parameter lines, list of bonds in in put molecule
    Output: Modified list of bond parameters
    Referenced By: GrabSmallMoleculeAMOEBAParameters
    Description:
    1. Iterate over list of bond parameter lines
    2. Grab poltype classes from parameter line
    3. For each poltype class, generate list of atom indices that have same type
    4. Now generate all combinations of indices that share same type
    5. Filter indices that arent consective in the molecule
    6. Grab the bond lengths for each bond in list of bonds with same types
    7. Average the bond lengths
    8. Replace bond length in parameter line with the average bond length 
    """
    newbondprms=[]
    # STEP 1
    for line in bondprms:
        linesplit=line.split()
        # STEP 2
        bondtypes=[int(linesplit[1]),int(linesplit[2])]
        bondindices=[]
        # STEP 3
        for prmtype in bondtypes:
            keylist=GrabKeysFromValue(poltype,poltype.idxtosymclass,prmtype)
            bondindices.append(keylist)
        # STEP 4
        allbonds = list(itertools.product(bondindices[0], bondindices[1]))
        allbonds=[x for x in allbonds if len(x) == len(set(x))]
        # STEP 5
        allbonds=FilterList(poltype,allbonds,bondlistbabel)
        tot=0
        # STEP 6
        for bond in allbonds:
            blen = optmol.GetBond(int(bond[0]),int(bond[1])).GetLength()
            tot+=blen
        if len(allbonds)==0:
            pass 
        else:
            # STEP 7
            avgbondlength=round(tot/len(allbonds),2)
            linesplit=re.split(r'(\s+)', line)   
            # STEP 8
            linesplit[8]=str(avgbondlength)
            line=''.join(linesplit)
        newbondprms.append(line)
    return newbondprms
      

def AddOptimizedAngleLengths(poltype,optmol,angleprms,anglelistbabel):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newangleprms=[]
    for line in angleprms:
        linesplit=line.split()
        angletypes=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
        angleindices=[]
        for prmtype in angletypes:
            keylist=GrabKeysFromValue(poltype,poltype.idxtosymclass,prmtype)
            angleindices.append(keylist)
        allangles = list(itertools.product(angleindices[0], angleindices[1],angleindices[2]))
        allangles=[x for x in allangles if len(x) == len(set(x))]
        allangles=FilterList(poltype,allangles,anglelistbabel)
        tot=0
        for angle in allangles:
            angle=[int(i) for i in angle]
            a = optmol.GetAtom(angle[0])
            b = optmol.GetAtom(angle[1])
            c = optmol.GetAtom(angle[2])
            anglelen = optmol.GetAngle(a,b,c)
            tot+=anglelen
        if len(allangles)==0:
            pass
        else:
            avganglelength=round(tot/len(allangles),2)
            linesplit=re.split(r'(\s+)', line)
            linesplit=linesplit[:11]
            linesplit.append('\n')
            transferredangle=float(linesplit[10])
            linesplit[10]=str(avganglelength)
            line=''.join(linesplit)
        line+='\n'
        newangleprms.append(line)
    return newangleprms
 


def GrabKeysFromValue(poltype,dic,thevalue):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    keylist=[]
    for key,value in dic.items():
        if value==thevalue:
            keylist.append(key)
    return keylist
          

def CheckIfAllTorsionsAreHydrogen(poltype,babelindices,mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    allhydrogentor=True
    atomobjects=[mol.GetAtom(i) for i in babelindices]
    a,b,c,d=atomobjects[:]
    aidx,bidx,cidx,didx=babelindices[:]
    aatomicnum=a.GetAtomicNum()
    datomicnum=d.GetAtomicNum()
    if aatomicnum!=1 or datomicnum!=1:
        allhydrogentor=False
    else:
        torlist=[]
        iteratomatom = openbabel.OBAtomAtomIter(b)
        for iaa in iteratomatom:
            iteratomatom2 = openbabel.OBAtomAtomIter(c)
            for iaa2 in iteratomatom2:
                ta = iaa.GetIdx()
                tb = bidx
                tc = cidx
                td = iaa2.GetIdx()
                if ((ta != tc and td != tb) and not (ta == aidx and td == didx)):
                    torlist.append([ta,tb,tc,td])
        for tor in torlist:
            atoms=[mol.GetAtom(i) for i in tor]
            aatomnum=atoms[0].GetAtomicNum()
            datomnum=atoms[3].GetAtomicNum()
            if aatomnum!=1 or datomnum!=1:
                allhydrogentor=False
    return allhydrogentor    

def CheckIfAllTorsionsAreHydrogenOneSide(poltype,babelindices,mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    allhydrogentoroneside=True
    atomobjects=[mol.GetAtom(i) for i in babelindices]
    a,b,c,d=atomobjects[:]
    aidx,bidx,cidx,didx=babelindices[:]
    aatomicnum=a.GetAtomicNum()
    datomicnum=d.GetAtomicNum()
    if aatomicnum!=1 and datomicnum!=1:
        allhydrogentoroneside=False
    else:
        torlist=[]
        
        iteratomatom = openbabel.OBAtomAtomIter(b)
        for iaa in iteratomatom:
            iteratomatom2 = openbabel.OBAtomAtomIter(c)
            for iaa2 in iteratomatom2:
                ta = iaa.GetIdx()
                tb = bidx
                tc = cidx
                td = iaa2.GetIdx()
                if ((ta != tc and td != tb) and not (ta == aidx and td == didx)):
                    torlist.append([ta,tb,tc,td])
        for tor in torlist:
            atoms=[mol.GetAtom(i) for i in tor]
            aatomnum=atoms[0].GetAtomicNum()
            datomnum=atoms[3].GetAtomicNum()
            if aatomnum!=1 and datomnum!=1:
                allhydrogentoroneside=False
    catomicnum=c.GetAtomicNum()
    batomicnum=b.GetAtomicNum()
    bhydcount=0
    iteratomatom = openbabel.OBAtomAtomIter(b)
    for iaa in iteratomatom:
        atomicnum=iaa.GetAtomicNum()
        if atomicnum==1:
            bhydcount+=1

    chydcount=0
    iteratomatom = openbabel.OBAtomAtomIter(c)
    for iaa in iteratomatom:
        atomicnum=iaa.GetAtomicNum()
        if atomicnum==1:
            chydcount+=1
    if batomicnum==6 and bhydcount==3:
        allhydrogentoroneside=False

    if catomicnum==6 and chydcount==3:
        allhydrogentoroneside=False



    return allhydrogentoroneside


def FindAllConsecutiveRotatableBonds(poltype,mol,listofbondsforprm):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    totalbondscollector=[]
    newrotbnds=[]
    for rotbnd in listofbondsforprm:
        b=rotbnd[0]+1
        c=rotbnd[1]+1
        batom=poltype.mol.GetAtom(b)
        catom=poltype.mol.GetAtom(c)
        bval=len([neighb for neighb in openbabel.OBAtomAtomIter(batom)])
        cval=len([neighb for neighb in openbabel.OBAtomAtomIter(catom)])
        if bval<2 or cval<2:
            continue
        newrotbnds.append(rotbnd)
    combs=list(combinations(newrotbnds,2)) 
    for comb in combs:
        firstbnd=comb[0]
        secondbnd=comb[1]
        total=firstbnd[:]+secondbnd[:]
        totalset=set(total)
        if len(totalset)==3:
            if firstbnd[-1]!=secondbnd[-1] and firstbnd[0]!=secondbnd[-1] and firstbnd[0]!=secondbnd[0]:
                secondbnd=secondbnd[::-1]
            elif firstbnd[0]==secondbnd[-1] and firstbnd[-1]!=secondbnd[-1] and firstbnd[0]!=secondbnd[0]:
                firstbnd=firstbnd[::-1]
            elif firstbnd[0]!=secondbnd[-1] and firstbnd[-1]!=secondbnd[-1] and firstbnd[0]==secondbnd[0]:
                firstbnd=firstbnd[::-1]
                secondbnd=secondbnd[::-1]


            catomindex=firstbnd[1]+1
            catom=poltype.mol.GetAtom(catomindex)
            isinring=catom.IsInRing()
            if isinring==True:
                continue
           
            totalbondscollector.append([firstbnd,secondbnd])
    return totalbondscollector

def FindAdjacentMissingTorsionsForTorTor(poltype,torsionsmissing,totalbondscollector,tortorsmissing):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for bndlist in totalbondscollector:
        first=bndlist[0]
        second=bndlist[1]
        foundfirst=CheckIfRotatableBondInMissingTorsions(poltype,first,torsionsmissing) 
        foundsecond=CheckIfRotatableBondInMissingTorsions(poltype,second,torsionsmissing) 
        foundtortormissing=CheckIfRotatableBondInMissingTorTors(poltype,first,second,tortorsmissing)
        if foundtortormissing==False:
            continue
        if (foundfirst==False and foundsecond==True and poltype.tortor==True):
            b=first[0]+1
            c=first[1]+1
            aatom,datom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(b),poltype.mol.GetAtom(c))
            a=aatom.GetIdx()
            d=datom.GetIdx()
            tor=[a-1,b-1,c-1,d-1]
            torsionsmissing.append(tuple(tor))
        elif (foundfirst==True and foundsecond==False and poltype.tortor==True):
            b=second[0]+1
            c=second[1]+1
            aatom,datom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(b),poltype.mol.GetAtom(c))
            a=aatom.GetIdx()
            d=datom.GetIdx()
            tor=[a-1,b-1,c-1,d-1]
            torsionsmissing.append(tuple(tor))
        elif (foundfirst==False and foundsecond==False and poltype.tortor==True):
            b=first[0]+1
            c=first[1]+1
            aatom,datom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(b),poltype.mol.GetAtom(c))
            a=aatom.GetIdx()
            d=datom.GetIdx()
            tor=[a-1,b-1,c-1,d-1]
            torsionsmissing.append(tuple(tor))
            b=second[0]+1
            c=second[1]+1
            aatom,datom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(b),poltype.mol.GetAtom(c))
            a=aatom.GetIdx()
            d=datom.GetIdx()
            tor=[a-1,b-1,c-1,d-1]
            torsionsmissing.append(tuple(tor))

    return torsionsmissing

def CheckIfRotatableBondInMissingTorsions(poltype,rotbnd,torsionsmissing):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    found=False
    for tor in torsionsmissing:
        bnd=[tor[1],tor[2]]
        if bnd==rotbnd or bnd[::-1]==rotbnd:
            found=True
            break
    return found

def CheckIfRotatableBondInMissingTorTors(poltype,firstbnd,secondbnd,tortormissing):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    foundtortormissing=False
    for tortor in tortormissing:
        first=[tortor[1],tortor[2]]
        second=[tortor[2],tortor[3]]
        if (first==firstbnd or first[::-1]==firstbnd) and (second==secondbnd or second[::-1]==secondbnd):
            foundtortormissing=True
        elif (first==secondbnd or first[::-1]==secondbnd) and (second==firstbnd or second[::-1]==firstbnd):
            foundtortormissing=True

    return foundtortormissing

def FindMissingTorTors(poltype,tortorindicestoextsmarts,tortorsmartsatomordertoparameters,rdkitmol,mol,indextoneighbidxs,totalbondscollector):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tortorsmissing=[]
    tortorsfound=[]
    for tortorindices,extsmarts in tortorindicestoextsmarts.items(): #only for ones that have match
        for tortorsmartsatomorder,parameters in tortorsmartsatomordertoparameters.items():
            smarts=tortorsmartsatomorder[0]
            if smarts==extsmarts:
                aidx,bidx,cidx,didx,eidx=tortorindices[:]
                firstneighborindexes=indextoneighbidxs[aidx]
                secondneighborindexes=indextoneighbidxs[bidx]
                thirdneighborindexes=indextoneighbidxs[cidx]
                fourthneighborindexes=indextoneighbidxs[didx]
                fifthneighborindexes=indextoneighbidxs[eidx]
                neighborindexes=firstneighborindexes+secondneighborindexes+thirdneighborindexes+fourthneighborindexes+fifthneighborindexes
                substructure = Chem.MolFromSmarts(smarts)
                matches=rdkitmol.GetSubstructMatches(substructure)
                matcharray=[]
                for match in matches:
                    for idx in match:
                        if idx not in matcharray:
                            matcharray.append(idx)
                check=CheckIfNeighborsExistInSMARTMatch(poltype,neighborindexes,matcharray)
                if check==False:
                    tortorsmissing.append(tortorindices)
                else:
                    tortorsfound.append(tortorindices)
    if poltype.onlyrotbndslist==None:
        poltype.onlyrotbndslist=[]
    for bndlist in totalbondscollector:
        first=bndlist[0]
        second=bndlist[1]
        babelfirst=[i+1 for i in first]
        babelsecond=[i+1 for i in second]
        #if not (babelfirst in poltype.onlyrotbndslist or babelfirst[::-1] in poltype.onlyrotbndslist):
        #    if (babelfirst in poltype.partialdoublebonds or babelfirst[::-1] in poltype.partialdoublebonds):
        #        continue
        #if not (babelsecond in poltype.onlyrotbndslist or babelsecond[::-1] in poltype.onlyrotbndslist):
        
        #    if (babelsecond in poltype.partialdoublebonds or babelsecond[::-1] in poltype.partialdoublebonds):
        #        continue    
        # allow partial double tor-tor until figure out non-adjacent tortor coupling
        b,c=first[:]
        d=second[0]
        aatom,dnewatom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(b+1),poltype.mol.GetAtom(c+1))
        bnewatom,eatom = torgen.find_tor_restraint_idx(poltype,poltype.mol,poltype.mol.GetAtom(c+1),poltype.mol.GetAtom(d+1))
        a=aatom.GetIdx()
        dnew=dnewatom.GetIdx()
        bnew=bnewatom.GetIdx()
        e=eatom.GetIdx()
        indices=[a-1,b,c,d,e-1]
        foundtortormissing=CheckIfRotatableBondInMissingTorTors(poltype,[bnew-1,c],[c,dnew-1],tortorsmissing)
        foundtortor=CheckIfRotatableBondInMissingTorTors(poltype,[bnew-1,c],[c,dnew-1],tortorsfound)
        if foundtortormissing==False and foundtortor==False:
            if len(poltype.onlyrottortorlist)==0:
                tortorsmissing.append(indices)
            else:
                ls=[b+1,c+1,d+1]
                if ls in poltype.onlyrottortorlist or ls[::-1] in poltype.onlyrottortorlist:
                    tortorsmissing.append(indices)

        else:
            ls=[b+1,c+1,d+1]
            if ls in poltype.onlyrottortorlist or ls[::-1] in poltype.onlyrottortorlist:
                tortorsmissing.append(indices)

    return tortorsmissing

def CheckIfRotatableBondsInOnlyRotBnds(poltype,first,second):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    inonlyrotbnds=True
    first=[i+1 for i in first]
    second=[i+1 for i in second]
    if (first in poltype.onlyrotbndslist or first[::-1] in poltype.onlyrotbndslist) and (second in poltype.onlyrotbndslist or second[::-1] in poltype.onlyrotbndslist):
        pass
    else:
        inonlyrotbnds=False

    return inonlyrotbnds 


def RingAtomicIndices(poltype,mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    sssr = mol.GetSSSR()
    atomindices=[]
    for ring in sssr:
        ringatomindices=GrabRingAtomIndices(poltype,mol,ring)
        atomindices.append(ringatomindices)        
    return atomindices

def GrabRingAtomIndices(poltype,mol,ring):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    ringatomindices=[]
    atomiter=openbabel.OBMolAtomIter(mol)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        if ring.IsInRing(atomidx)==True:
            ringatomindices.append(atomidx)
    return ringatomindices

def GrabRingAtomIndicesFromInputIndex(poltype,atomindex,atomindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for ring in atomindices:
        if atomindex in ring:
            return ring

    ring=None
    return ring

def GrabIndicesInRing(poltype,babelindices,ring):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    ringtorindices=[]
    for index in ring:
        if index in babelindices:
            ringtorindices.append(index)

    return ringtorindices


def FindMissingTorsions(poltype,torsionindicestoparametersmartsenv,rdkitmol,mol,indextoneighbidxs):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    torsionsmissing=[]
    poormatchingaromatictorsions=[]
    poormatchingpartialaromatictorsions=[]
    torsionstozerooutduetocolinear=[]
    for torsionindices,smartsenv in torsionindicestoparametersmartsenv.items():
        if poltype.dontdotor==True:
            continue
        aidx,bidx,cidx,didx=torsionindices[:]
        babelindices=[i+1 for i in torsionindices]

        allhydrogentor=CheckIfAllTorsionsAreHydrogen(poltype,babelindices,mol)
        allhydrogentoroneside=CheckIfAllTorsionsAreHydrogenOneSide(poltype,babelindices,mol)
        atoma=rdkitmol.GetAtomWithIdx(aidx)
        atomd=rdkitmol.GetAtomWithIdx(didx)
        atomicnumatoma=atoma.GetAtomicNum()
        atomicnumatomd=atomd.GetAtomicNum()
        abidx,bbidx,cbidx,dbidx=babelindices[:]
        
        bond=mol.GetBond(bbidx,cbidx)
        bondorder=bond.GetBondOrder()
        babelatoms=[mol.GetAtom(i) for i in babelindices]
        aatom,batom,catom,datom=babelatoms[:]

        middlebond=mol.GetBond(babelindices[1],babelindices[2])
        ringbond=middlebond.IsInRing()
        firstangle=mol.GetAngle(aatom,batom,catom)
        secondangle=mol.GetAngle(batom,catom,datom)
        if firstangle<0:
            firstangle=firstangle+360
        if secondangle<0:
            secondangle=secondangle+360
        angletol=3.5
        if np.abs(180-firstangle)<=angletol or np.abs(180-secondangle)<=angletol:
            torsionstozerooutduetocolinear.append(torsionindices)
            continue
         
        
        atomvals=[len([neighb for neighb in openbabel.OBAtomAtomIter(a)]) for a in babelatoms]
        atomnums=[a.GetAtomicNum() for a in babelatoms]
        batomnum=atomnums[1]
        catomnum=atomnums[2]
           
        ringbools=[a.IsInRing() for a in babelatoms]
        arobools=[a.IsAromatic() for a in babelatoms]
        bnd=[babelindices[1],babelindices[2]]
         
        ringb=ringbools[1]
        ringc=ringbools[2]
        aroa=arobools[0]
        arob=arobools[1]
        aroc=arobools[2]
        arod=arobools[3]
        hybs=[a.GetHyb() for a in babelatoms]
        hybb=hybs[1]
        hybc=hybs[2]
        firstneighborindexes=indextoneighbidxs[aidx]
        secondneighborindexes=indextoneighbidxs[bidx]
        thirdneighborindexes=indextoneighbidxs[cidx]
        fourthneighborindexes=indextoneighbidxs[didx]
        neighborindexes=firstneighborindexes+secondneighborindexes+thirdneighborindexes+fourthneighborindexes
        smarts=smartsenv[0]
        substructure = Chem.MolFromSmarts(smarts)
        matches=rdkitmol.GetSubstructMatches(substructure)
        for match in matches:
            allin=True
            for idx in torsionindices:
                if idx not in match:
                    allin=False
            if allin==True:
                matcharray=match
                break

        check=CheckIfNeighborsExistInSMARTMatch(poltype,neighborindexes,matcharray)
        if len(poltype.onlyrotbndslist)!=0:
            if [bbidx,cbidx] in poltype.onlyrotbndslist or [cbidx,bbidx] in poltype.onlyrotbndslist:
                check=False 
        if '~' in smarts or '*' in smarts:
            check=False
        if [bbidx,cbidx] in poltype.dontrotbndslist or [cbidx,bbidx] in poltype.dontrotbndslist:
            check=True
        if ringbond==True:
            atomindices=RingAtomicIndices(poltype,mol)
            therings=torgen.GrabAllRingsContainingMostIndices(poltype,atomindices,babelindices,2)
            if (len(therings)>0) and poltype.dontfrag==False:
                if len(therings[0])>7: # special case where whole molecule is a ring then dont consider ring bond
                    if hybs[1]!=2 and hybs[2]!=2:
                        ringbond=False

            ring=GrabRingAtomIndicesFromInputIndex(poltype,babelindices[1],atomindices)
            ringtorindices=GrabIndicesInRing(poltype,babelindices,ring)
        if (bnd in poltype.partialdoublebonds or bnd[::-1] in poltype.partialdoublebonds) and poltype.rotalltors==False and ([bbidx,cbidx] not in poltype.onlyrotbndslist and [cbidx,bbidx] not in poltype.onlyrotbndslist) and check==True:
            continue

        if check==False:
            if ringbond==True:
                if (2 not in hybs): # non-aromatic torsion want parameters for 
                    if poltype.transferanyhydrogentor==True and (atomicnumatoma==1 or atomicnumatomd==1) and (allhydrogentor==False and allhydrogentoroneside==False): # then here transfer torsion because can pick up most QM-MM on heavy atoms, less parameters to fit
                        poormatchingpartialaromatictorsions.append(torsionindices)
                    else: # if dont have heavy atoms on either side then just fit the hydrogen torsion
                        if poltype.nonaroringtor1Dscan==True or poltype.refinenonaroringtors==True: 
                            found=False
                            if len(poltype.onlyrotbndslist)!=0:
                                if [bbidx,cbidx] in poltype.onlyrotbndslist or [cbidx,bbidx] in poltype.onlyrotbndslist:
                                    found=True
                            else:
                                found=True
                            if found==True:
                                if len(ring)>3:
                                    if torsionindices not in torsionsmissing and poltype.dontfrag==False: # make sure fragmenter is on (wont work for < 25 atoms by default)
                                        torsionsmissing.append(torsionindices)
                            else:
                                poormatchingpartialaromatictorsions.append(torsionindices)
                        else:
                            poormatchingpartialaromatictorsions.append(torsionindices)

                        
                elif hybs[1]==2 and hybs[2]==2:
                    if torsionindices not in poormatchingaromatictorsions:
                        poormatchingaromatictorsions.append(torsionindices)
                elif (hybs[1]==2 and hybs[2]!=2) or (hybs[1]!=2 and hybs[2]==2):
                    if torsionindices not in poormatchingpartialaromatictorsions:
                        poormatchingpartialaromatictorsions.append(torsionindices)
                elif (hybs[1]!=2 and hybs[2]!=2) and (hybs[0]==2 or hybs[3]==2):
                    if torsionindices not in poormatchingpartialaromatictorsions:
                        poormatchingpartialaromatictorsions.append(torsionindices)

            else:
                if poltype.transferanyhydrogentor==True and poltype.rotalltors==False and (atomicnumatoma==1 or atomicnumatomd==1) and (allhydrogentor==False and allhydrogentoroneside==False): # then here transfer torsion because can pick up most QM-MM on heavy atoms, less parameters to fit
                    if hybs[1]==2 and hybs[2]==2:
                        if torsionindices not in poormatchingaromatictorsions:
                            poormatchingaromatictorsions.append(torsionindices)
                    elif (hybs[1]==2 and hybs[2]!=2) or (hybs[1]!=2 and hybs[2]==2):
                        if torsionindices not in poormatchingpartialaromatictorsions:
                            poormatchingpartialaromatictorsions.append(torsionindices)
                    elif (hybs[1]!=2 and hybs[2]!=2) and (hybs[0]==2 or hybs[3]==2):
                        if torsionindices not in poormatchingpartialaromatictorsions:
                            poormatchingpartialaromatictorsions.append(torsionindices)
                    else:
                        if torsionindices not in poormatchingpartialaromatictorsions:
                            poormatchingpartialaromatictorsions.append(torsionindices)


                else:
                    if len(poltype.onlyrotbndslist)!=0:
                        if [bbidx,cbidx] in poltype.onlyrotbndslist or [cbidx,bbidx] in poltype.onlyrotbndslist:
                            if bondorder!=1: 
                                continue

                            if torsionindices not in torsionsmissing:
                                torsionsmissing.append(torsionindices)
                    else:
                        if bondorder!=1:
                            continue
                        if torsionindices not in torsionsmissing:
                            torsionsmissing.append(torsionindices)

                    continue

        else:
            if poltype.rotalltors==True and ringbond==False:
                if bondorder<2:
                    if torsionindices not in torsionsmissing:
                        torsionsmissing.append(torsionindices)
            if len(poltype.onlyrotbndslist)!=0:
                if [bbidx,cbidx] in poltype.onlyrotbndslist or [cbidx,bbidx] in poltype.onlyrotbndslist:
                    if bondorder<2:
                        if torsionindices not in torsionsmissing:
                            torsionsmissing.append(torsionindices)
    return torsionsmissing,poormatchingaromatictorsions,poormatchingpartialaromatictorsions,torsionstozerooutduetocolinear


def GrabAllRingsContainingIndices(poltype,atomindices,babelindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    rings=[]
    for ring in atomindices:
        allin=True
        for atomindex in babelindices:
            if atomindex not in ring:
                allin=False
        if allin==True:
            rings.append(ring)
    return rings

def FindAllNeighborIndexes(poltype,rdkitmol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indextoneighbidxs={}
    for atm in rdkitmol.GetAtoms():
        atmidx=atm.GetIdx()
        if atmidx not in indextoneighbidxs.keys():
            indextoneighbidxs[atmidx]=[]
        for neighbatm in atm.GetNeighbors():
            neighbatmidx=neighbatm.GetIdx()
            if neighbatmidx not in indextoneighbidxs[atmidx]:
                indextoneighbidxs[atmidx].append(neighbatmidx)
            

    return indextoneighbidxs

def CheckIfNeighborsExistInSMARTMatch(poltype,neighborindexes,smartsindexes):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    check=True
    for idx in neighborindexes:
        if idx not in smartsindexes:
            check=False
    return check

def ZeroOutMissingStrbnd(poltype,anglemissingtinkerclassestopoltypeclasses,strbndprms):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newstrbndprms=[]
    for line in strbndprms:
        linesplit=line.split()
        classes=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3])])
        for sublist in anglemissingtinkerclassestopoltypeclasses.values():
            if classes in sublist or classes[::-1] in sublist: 
                newlinesplit=re.split(r'(\s+)', line)
                newlinesplit[8]='0'
                newlinesplit[10]='0'
                line=''.join(newlinesplit)
        newstrbndprms.append(line)

    return newstrbndprms

def AssignAngleGuessParameters(poltype,anglemissingtinkerclassestopoltypeclasses,angleprms,indextoneighbidxs):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newangleprms=[]

    for line in angleprms:
        linesplit=line.split()
        classes=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3])])
        for tinkerclasses,sublist in anglemissingtinkerclassestopoltypeclasses.items():
            found=False
            if classes in sublist:
                found=True
            elif classes[::-1] in sublist: 
                found=True
            if found==True:
                indices=[]
                for poltypeclass in classes:
                    keylist=GrabKeysFromValue(poltype,poltype.idxtosymclass,poltypeclass)
                    indices.append(keylist)
                exampleindices=None
                combs = list(itertools.product(*indices))
                for comb in combs:
                    poscomb=[np.abs(i) for i in comb]
                    poscomb=[k-1 for k in poscomb]
                    checkconsec=CheckIfAtomsConnected(poltype,poscomb,indextoneighbidxs)
                    if checkconsec==True:
                        exampleindices=[k-1 for k in comb]
                        break
                exampleindices=[int(i) for i in exampleindices]
                atoms=[poltype.rdkitmol.GetAtomWithIdx(k) for k in exampleindices]
                atomicnums=[a.GetAtomicNum() for a in atoms]
                atomicval=[a.GetExplicitValence() for a in atoms]

                angleguess=AngleGuess(poltype,atomicnums[0],atomicnums[1],atomicnums[2],atomicval[0],atomicval[1],atomicval[2])
                newlinesplit=re.split(r'(\s+)', line)
                newlinesplit[8]=str(angleguess)
                line=''.join(newlinesplit)
        newangleprms.append(line)

    return newangleprms


def CheckIfAtomsConnected(poltype,poscomb,endindextoneighbs):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    checkconsec=True
    if len(poscomb)>1:
        for i in range(len(poscomb)-1):
            index=poscomb[i]
            nextindex=poscomb[i+1]
            indexneighbs=endindextoneighbs[index]
            if nextindex not in indexneighbs:
                checkconsec=False



    return checkconsec

def AssignBondGuessParameters(poltype,bondmissingtinkerclassestopoltypeclasses,bondprms,indextoneighbidxs):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newbondprms=[]
    for line in bondprms:
        linesplit=line.split()
        classes=tuple([int(linesplit[1]),int(linesplit[2])])
        for tinkerclasses,sublist in bondmissingtinkerclassestopoltypeclasses.items():
            found=False
            if classes in sublist:
                found=True
            elif classes[::-1] in sublist: 
                found=True
            if found==True:
                indices=[]
                for poltypeclass in classes:
                    keylist=GrabKeysFromValue(poltype,poltype.idxtosymclass,poltypeclass)
                    indices.append(keylist)
                exampleindices=None
                combs = list(itertools.product(*indices))
                for comb in combs:
                    poscomb=[np.abs(i) for i in comb]
                    poscomb=[k-1 for k in poscomb]
                    checkconsec=CheckIfAtomsConnected(poltype,poscomb,indextoneighbidxs)
                    if checkconsec==True:
                        exampleindices=[k-1 for k in comb]
                        break 
                exampleindices=[int(i) for i in exampleindices]
                atoms=[poltype.rdkitmol.GetAtomWithIdx(k) for k in exampleindices]
                atomicnums=[a.GetAtomicNum() for a in atoms]
                atomicval=[a.GetExplicitValence() for a in atoms]

                bondguess=BondGuess(poltype,atomicnums[0],atomicnums[1],atomicval[0],atomicval[1])
                newlinesplit=re.split(r'(\s+)', line)
                newlinesplit[6]=str(bondguess)
                line=''.join(newlinesplit)
        newbondprms.append(line)

    return newbondprms




def ZeroOutMissingTorsions(poltype,torsionsmissingtinkerclassestopoltypeclasses,torsionprms):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newtorsionprms=[]
    for line in torsionprms:
        linesplit=line.split()
        classes=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3]),int(linesplit[4])])
        for sublist in torsionsmissingtinkerclassestopoltypeclasses.values():
            if classes in sublist or classes[::-1] in sublist: 
                newlinesplit=re.split(r'(\s+)', line)
                splitafter=newlinesplit[10:]
                splitafter[0]='0'
                splitafter[6]='0'
                splitafter[12]='0'
                newlinesplit=newlinesplit[:10]+splitafter
                line=''.join(newlinesplit)
        newtorsionprms.append(line)

    return newtorsionprms


def GrabKeysFromValue(poltype,dic,thevalue):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    keylist=[]
    for key,value in dic.items():
        if value==thevalue:
            keylist.append(key)
    return keylist


def DefaultAromaticMissingTorsions(poltype,arotorsionsmissingtinkerclassestopoltypeclasses,partialarotorsionsmissingtinkerclassestopoltypeclasses,torsionprms,mol): # transfer bezene aromatic torsion paramerers from amoeba09
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newtorsionprms=[]
    poorarohydcounttoprms={0:[-.67,6.287,0],1:[.55,6.187,-.55],2:[0,6.355,0]}
    poorpartialarohydcounttoprms={0:[.854,-.374,.108],1:[0,0,.108],2:[0,0,.299]}
    poorarohydcounttodescrips={0:["Benzene C","Benzene C","Benzene C","Benzene C"],1:["Benzene HC","Benzene C","Benzene C","Benzene C"],2:["Benzene HC","Benzene C","Benzene C","Benzene HC"]}
    poorpartialarohydcounttodescrips={0:["Alkane -CH2-","Alkane -CH2-","Alkane -CH2-","Alkane -CH2-"],1:["Alkane -H2C-","Alkane -CH2-","Alkane -CH2-","Alkane -CH2-"],2:["Alkane -H2C-","Alkane -CH2-","Alkane -CH2-","Alkane -H2C-"]}
    ringextra='# Ring bond detected for missing torsion'+'\n'
    doubleextra='# Double bond detected for missing torsion'+'\n'
    singleextra='# Single bond detected for missing torsion'+'\n'
    hydextra='# Transferring hydrogen torsion to reduce fit parameters'+'\n'
    arotorsionlinetodescrips={} 
    hydclasses=[]
    for atom in poltype.rdkitmol.GetAtoms():
        atomicnum=atom.GetAtomicNum()
        babelindex=atom.GetIdx()+1
        symclass=poltype.idxtosymclass[babelindex]
        if atomicnum==1:
            hydclasses.append(symclass)
    
    for line in torsionprms:
        linesplit=line.split()
        a=int(linesplit[1])
        b=int(linesplit[2])
        c=int(linesplit[3])
        d=int(linesplit[4])
        allcombs=[]
        ls=[a,b,c,d]
        for typenum in ls: 
            indices=GrabKeysFromValue(poltype,poltype.idxtosymclass,typenum)
            allcombs.append(indices)
        combs=list(itertools.product(*allcombs))
        for comb in combs:
            aindex=comb[0]
            bindex=comb[1]
            cindex=comb[2]
            dindex=comb[3]
            firstbond=mol.GetBond(aindex,bindex)
            midbond=mol.GetBond(bindex,cindex)
            lastbond=mol.GetBond(cindex,dindex)
            if midbond!=None and firstbond!=None and lastbond!=None:
                break 
        babelindices=list(comb)
        babelatoms=[mol.GetAtom(i) for i in babelindices]
        hybs=[a.GetHyb() for a in babelatoms]
        ringbond=midbond.IsInRing()
        if ringbond==True:
            atomindices=RingAtomicIndices(poltype,mol)
            therings=torgen.GrabAllRingsContainingMostIndices(poltype,atomindices,babelindices,2)
            if (len(therings)>0) and poltype.dontfrag==False:
                if len(therings[0])>7: # special case where whole molecule is a ring then dont consider ring bond
                    if hybs[1]!=2 and hybs[2]!=2:
                        ringbond=False

        bondorder=midbond.GetBondOrder()      
        classes=tuple([a,b,c,d])
        hydcount=0
        if a in hydclasses:
            hydcount+=1
        if d in hydclasses:
            hydcount+=1
        found=False
        for sublist in arotorsionsmissingtinkerclassestopoltypeclasses.values():
            if classes in sublist or classes[::-1] in sublist: 
                found=True
                prms=poorarohydcounttoprms[hydcount]
                descrips=poorarohydcounttodescrips[hydcount]
                break
        for sublist in partialarotorsionsmissingtinkerclassestopoltypeclasses.values():
            if classes in sublist or classes[::-1] in sublist: 
                found=True
                prms=poorpartialarohydcounttoprms[hydcount]
                descrips=poorpartialarohydcounttodescrips[hydcount]
                break

        if found==True:
            newlinesplit=re.split(r'(\s+)', line)
            splitafter=newlinesplit[10:]
            splitafter[0]=str(prms[0])
            splitafter[6]=str(prms[1])
            splitafter[12]=str(prms[2])
            newlinesplit=newlinesplit[:10]+splitafter
            line=''.join(newlinesplit)
            if ringbond==True:
                extra=ringextra
            elif ringbond==False and bondorder!=1:
                extra=doubleextra
            else:
                extra=singleextra
            if hydcount>0 and ringbond==False:
                extra+=hydextra
            arotorsionlinetodescrips[line]=[extra,descrips]
        newtorsionprms.append(line)

    return newtorsionprms,arotorsionlinetodescrips


def GrabTorsionParameterCoefficients(poltype,torsionprms):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    torsionkeystringtoparameters={}
    for line in torsionprms:
        linesplit=line.split() 
        key = '%s %s %s %s' % (linesplit[1], linesplit[2], linesplit[3], linesplit[4])
        parameters=[float(linesplit[5]),float(linesplit[8]),float(linesplit[11])]
        allzero=True
        for prm in parameters:
            if prm!=0:
                allzero=False
        if allzero==False:
            torsionkeystringtoparameters[key]=parameters
       
    return torsionkeystringtoparameters

def PruneDictionary(poltype,keysubset,dic):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newdic={}
    removeditems=[]
    for key in keysubset:
        if key in dic.keys():
            value=dic[key]
            newdic[key]=value

        elif key[::-1] in dic.keys():
            value=dic[key[::-1]]
            newdic[key]=value
        else:
            removeditems.append(key)

    return newdic,removeditems


def TinkerClassesToPoltypeClasses(poltype,indicestotinkerclasses,formatorder=True):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tinkerclassestopoltypeclasses={}
    poltypeclassesalreadyassigned=[]
    for indices,tinkerclasses in indicestotinkerclasses.items():
        babelindices=[i+1 for i in indices]
        poltypeclasses=[poltype.idxtosymclass[i] for i in babelindices]
        revtinkerclasses=tinkerclasses[::-1]
        if poltypeclasses in poltypeclassesalreadyassigned:
            continue
        if formatorder==True:
            if len(indices)>1:
                first=int(tinkerclasses[0])
                second=int(tinkerclasses[-1])
                if first>second:
                    pass
                elif first<second:
                    tinkerclasses=tinkerclasses[::-1]
                    poltypeclasses=poltypeclasses[::-1]
                else:
                    if len(indices)==4:
                        first=int(tinkerclasses[0])
                        second=int(tinkerclasses[1])
                        third=int(tinkerclasses[2])
                        fourth=int(tinkerclasses[3])
                        firstsum=first+second
                        secondsum=third+fourth
                        if firstsum>secondsum:
                            pass
                        elif firstsum<secondsum:
                            tinkerclasses=tinkerclasses[::-1]
                            poltypeclasses=poltypeclasses[::-1]

        poltypeclassesalreadyassigned.append(poltypeclasses)                
        if tuple(tinkerclasses) not in tinkerclassestopoltypeclasses.keys():
            tinkerclassestopoltypeclasses[tuple(tinkerclasses)]=[]
        if tuple(poltypeclasses) not in tinkerclassestopoltypeclasses[tuple(tinkerclasses)]: 
            tinkerclassestopoltypeclasses[tuple(tinkerclasses)].append(tuple(poltypeclasses))

    return tinkerclassestopoltypeclasses


def ConvertIndicesDictionaryToPoltypeClasses(poltype,indicestovalue,indicestotinkerclasses,tinkerclassestopoltypeclasses):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltypeclassestovalue={}
    for indices,value in indicestovalue.items():
        if indices in indicestotinkerclasses.keys():
            tinkerclasses=tuple(indicestotinkerclasses[indices])
            if tinkerclasses in tinkerclassestopoltypeclasses.keys():
                poltypeclasses=tuple(tinkerclassestopoltypeclasses[tinkerclasses])
            elif tinkerclasses[::-1] in tinkerclassestopoltypeclasses.keys():
                poltypeclasses=tuple(tinkerclassestopoltypeclasses[tinkerclasses[::-1]])
            else:
                continue
            babelindices=[i+1 for i in indices]
            symclasses=tuple([poltype.idxtosymclass[i] for i in babelindices])
            for ls in poltypeclasses:
                if ls==symclasses or ls==symclasses[::-1]:
                    ls=tuple([ls])
                    if ls not in poltypeclassestovalue.keys():
                        poltypeclassestovalue[ls]=[]
                        if value not in poltypeclassestovalue[ls]:
                            poltypeclassestovalue[ls].append(value)
    return poltypeclassestovalue 

def CheckIfStringInStringList(poltype,string,stringlist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    found=False
    for e in stringlist:
        if e==string:
            found=True
    return found

def GrabTypesFromPrmLine(poltype,ls):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    typenums=[]
    for e in ls:
        isint=RepresentsInt(e)  
        if isint==True:
            typenums.append(e) 
    if len(typenums)>=4:
        if typenums[-2]=='0' and typenums[-1]=='0':
            typenums=typenums[:2]
    if len(typenums)>4:
        typenums=typenums[:4]     
    return typenums


def SearchForPoltypeClasses(poltype,prmline,poltypeclasseslist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    listofpoltypeclasses=None
    poltypeclasses=None
    allin=None
    prmlinesplit=prmline.split()
    typenums=GrabTypesFromPrmLine(poltype,prmlinesplit)
    for listofpoltypeclasses in poltypeclasseslist:
        for poltypeclasses in listofpoltypeclasses:
            allin=True
            if int(typenums[0])!=poltypeclasses[0]: # then try flipping
                poltypeclasses=poltypeclasses[::-1]
            for i in range(len(typenums)):
                poltypeclass=int(typenums[i])
                otherpoltypeclass=poltypeclasses[i]
                if poltypeclass!=otherpoltypeclass:
                    allin=False 
            if allin==True:
                return listofpoltypeclasses,poltypeclasses
    if allin!=None:
        if allin==False:
            listofpoltypeclasses=None
            poltypeclasses=None

    return listofpoltypeclasses,poltypeclasses


def MapParameterLineToTransferInfo(poltype,prms,poltypeclassestoparametersmartsatomorders,poltypeclassestosmartsatomorders,poltypeclassestoelementtinkerdescrips,poltypeclassestosmartsatomordersext,newpoltypeclassestocomments,newpoltypeclassestosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses,defaultvalues=None,keyword=None):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    prmstotransferinfo={}
    for line in prms:
        warn=False
        if keyword!=None:
            if keyword not in line:
                transferinfo='# blank'
            else:
                poltypeclasses,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,poltypeclassestosmartsatomordersext.keys())
                smartsatomorders=poltypeclassestosmartsatomordersext[poltypeclasses]
                transferinfoline='#'+' '+'matching SMARTS from molecule '+' '+str(smartsatomorders)+' from external database'+'\n'

        else:
            if line in arotorsionlinetodescrips.keys():
                tup=arotorsionlinetodescrips[line]
                extra=tup[0]
                descrips=tup[1]
                transferinfoline=extra+'# Transferring from '+str(descrips)+'\n' 
            else:
                poltypeclasses,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,poltypeclassestoparametersmartsatomorders.keys())
                if poltypeclasses==None:
                         
                    poltypeclasses,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,poltypeclassestosmartsatomordersext.keys())

                    if poltypeclasses==None:
                        poltypeclasses,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,newpoltypeclassestocomments.keys())

                        comments=newpoltypeclassestocomments[poltypeclasses]
                        smartslist=newpoltypeclassestosmartslist[poltypeclasses]
                        commentstring=' '.join(comments)
                        smartsliststring=' '.join(smartslist)
                        transferinfoline='# amoeba21'+' '+'comments='+commentstring+' '+'SMARTS match = '+smartsliststring+'\n'
                    else:
                        smartsatomorders=poltypeclassestosmartsatomordersext[poltypeclasses]
                        transferinfoline='#'+' '+'matching SMARTS from molecule '+' '+str(smartsatomorders)+' from external database'+'\n'

                else:
                    parametersmartsatomorders=poltypeclassestoparametersmartsatomorders[poltypeclasses]
                    smartsatomorders=poltypeclassestosmartsatomorders[poltypeclasses]
                    elementtinkerdescrips=poltypeclassestoelementtinkerdescrips[poltypeclasses]
                    transferinfoline='# amoeba09'+' '+'matching SMARTS from molecule '+' '+str(smartsatomorders)+' '+'to SMARTS from parameter file'+' '+str(parametersmartsatomorders)+' '+'with tinker type descriptions '+str(elementtinkerdescrips)+'\n'
                    if defaultvalues!=None:
                        for value in defaultvalues:
                            if str(value) in line:
                                warn=True
                    if warn==True:
                        transferinfoline+='# '+'WARNING DEFAULT MM3 OPBEND VALUES USED '+'\n'
                    smarts=smartsatomorders[0]
                    if '~' in smarts or '*' in smarts:
                        transferinfoline+='# '+'WARNING WILDCARDS USED IN SMARTS PARAMETER MATCHING'+'\n'
                        warn=True
        showtransferinfo=True
        extraline=''
        if 'vdw' in line:
            linesplit=line.split()
            vdwtype=int(linesplit[1])
            if vdwtype in missingvdwtypes:
                warn=True
                extraline+='# Missing vdw parameters'+'\n'

        if 'angle' in line:
            linesplit=line.split()
            a=int(linesplit[1])
            b=int(linesplit[2])
            c=int(linesplit[3])
            ls=[a,b,c]
            if ls in missinganglepoltypeclasses or ls[::-1] in missinganglepoltypeclasses:
                extraline+='# Missing angle parameters, assigning default parameters via element and valence'+'\n'
                warn=True
                showtransferinfo=False


        if 'strbnd' in line:
            linesplit=line.split()
            a=int(linesplit[1])
            b=int(linesplit[2])
            c=int(linesplit[3])
            ls=[a,b,c]
            if ls in missinganglepoltypeclasses or ls[::-1] in missinganglepoltypeclasses:
                extraline+='# Missing strbnd parameters, zeroing out parameters'+'\n'
                warn=True
                showtransferinfo=False

        if 'bond' in line:
            linesplit=line.split()
            a=int(linesplit[1])
            b=int(linesplit[2])
            ls=[a,b]
            if ls in missingbondpoltypeclasses or ls[::-1] in missingbondpoltypeclasses:
                extraline+='# Missing bond parameters, assigning default parameters via element and valence'+'\n'
                warn=True
                showtransferinfo=False


        if 'torsion' in line:
            linesplit=line.split()
            a=int(linesplit[1])
            b=int(linesplit[2])
            c=int(linesplit[3])
            d=int(linesplit[4])
            ls=[a,b,c,d]
            if ls in torsionsmissing or ls[::-1] in torsionsmissing:
                extraline+='# Missing torsion parameters, will attempt to fit parameters'+'\n'
                showtransferinfo=True


        transferinfoline+=extraline
        if showtransferinfo==False:
            transferinfoline=extraline
        if warn==True:
            poltype.WriteToLog(transferinfoline)
            warnings.warn(transferinfoline)
        prmstotransferinfo[line]=transferinfoline

    return prmstotransferinfo

def FindMaximumCommonSubstructures(poltype,parametersmartslist,rdkitmol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    parametersmartstomaxcommonsubstructure={}
    maxatomsize=0
    for parametersmarts in parametersmartslist:
        mols = [rdkitmol,Chem.MolFromSmarts(parametersmarts)]
        res=rdFMCS.FindMCS(mols)
        atomnum=res.numAtoms
        smartsmcs=res.smartsString
        if atomnum>0:
            if atomnum>maxatomsize:
                maxatomsize=atomnum
            parametersmartstomaxcommonsubstructure[parametersmarts]=smartsmcs
    return parametersmartstomaxcommonsubstructure,maxatomsize 

def GrabPlanarBonds(poltype,listofbondsforprm,mol): # used for checking missing opbend parameters later
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    planarbonds=[]
    for bond in listofbondsforprm:
        revbond=bond[::-1]
        totalbonds=[bond,revbond]
        for lsidx in range(len(totalbonds)):
            bond=totalbonds[lsidx]
            a = mol.GetAtom(bond[0]+1)
            b = mol.GetAtom(bond[1]+1)
            ainring=a.IsInRing()
            binring=b.IsInRing()
            aisaromatic=a.IsAromatic()
            bisaromatic=b.IsAromatic()
            if ainring==True and binring==True:
                if aisaromatic==True or bisaromatic==True:
                    if bond not in planarbonds:
                        planarbonds.append(bond)
                else:
                    if b.GetHyb()==2 and len(list(openbabel.OBAtomAtomIter(b)))==3:
                        if bond not in planarbonds:
                            planarbonds.append(bond)

            else:
                if b.GetHyb()==2 and len(list(openbabel.OBAtomAtomIter(b)))==3:
                    if bond not in planarbonds:
                        planarbonds.append(bond)
    return planarbonds 


def FindPotentialMissingParameterTypes(poltype,prms,tinkerclassestopoltypeclasses):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltypeclasseslist=list(tinkerclassestopoltypeclasses.values())
    newpoltypeclasseslist=[]
    for poltypeclassesls in poltypeclasseslist:
        newpoltypeclassesls=[]
        for poltypeclasses in poltypeclassesls:
            newpoltypeclassesls.append(poltypeclasses)
        newpoltypeclasseslist.append(newpoltypeclassesls)
    missingprms=[]
    foundprms=[]
    for line in prms:
        poltypeclassesls,subpoltypeclasses=SearchForPoltypeClasses(poltype,line,newpoltypeclasseslist)
        linesplit=line.split()
        prms=[int(linesplit[1]),int(linesplit[2])]
        if poltypeclassesls!=None:
            for poltypeclasses in poltypeclassesls:
                poltypeclasses=tuple([int(i) for i in poltypeclasses])
                if poltypeclasses[0]==prms[0] and poltypeclasses[1]==prms[1]:
                    foundprms.append(poltypeclasses) 
    for poltypeclassesls in newpoltypeclasseslist:
        for poltypeclasses in poltypeclassesls:
            if poltypeclasses not in foundprms and poltypeclasses not in missingprms:
                missingprms.append(poltypeclasses)
 
    return missingprms 


def ConvertPoltypeClassesToIndices(poltype,missingprmtypes):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    missingprmindices=[]
    for poltypeclasses in missingprmtypes:
        indices=[]
        for poltypeclass in poltypeclasses:
            keylist=GrabKeysFromValue(poltype,poltype.idxtosymclass,poltypeclass)
            keylist=[i-1 for i in keylist]
            indices.append(keylist)
        missingprmindices.append(indices)
    finallist=[]
    for indices in missingprmindices:
        combs = list(itertools.product(*indices))
        for comb in combs:
            finallist.append(comb)
    return finallist

def FilterIndices(poltype,potentialmissingprmindices,indices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    missingprmindices=[]
    for prmindices in potentialmissingprmindices:
        if list(prmindices) in indices:
            missingprmindices.append(prmindices)
    return missingprmindices
        

def GrabTypeAndDescriptions(poltype,prmfile):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    elementmm3descriptomm3type={}
    temp=open(prmfile,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'atom' in line and 'piatom' not in line:
            linesplit=line.split()    
            newlinesplit=linesplit[1:-3]
            mm3type=newlinesplit[0]
            element=newlinesplit[1]
            mm3descrip=' '.join(newlinesplit[2:])
            ls=[element,mm3descrip]
            elementmm3descriptomm3type[tuple(ls)]=mm3type
    return elementmm3descriptomm3type 


def DefaultOPBendParameters(poltype,missingopbendprmindices,mol,opbendbondindicestotrigonalcenterbools):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newopbendprms=[]
    defaultvalues=[]
    for opbendprmindices in missingopbendprmindices:
        babelindices=[i+1 for i in opbendprmindices]
        poltypeclasses=[poltype.idxtosymclass[i] for i in babelindices]
        boolarray=opbendbondindicestotrigonalcenterbools[opbendprmindices]
        if boolarray[1]==True:
            atomindex=int(babelindices[1])
            atom=mol.GetAtom(atomindex)
            neighbs=list(openbabel.OBAtomAtomIter(atom))
            if len(neighbs)==3:
                atomnum=atom.GetAtomicNum()
                if atomnum==6:
                    match=False
                    for natom in neighbs:
                        natomicnum=natom.GetAtomicNum()
                        if natomicnum==8:
                            bnd=mol.GetBond(atom,natom)
                            BO=bnd.GetBondOrder()
                            if BO==2:
                                match=True
                    if match==True:
                        opbendvalue=round(1*71.94,2)
                    else:
                        opbendvalue=round(.2*71.94,2)
                elif atomnum==7:
                    count=0
                    for natom in neighbs:
                        natomicnum=natom.GetAtomicNum()
                        if natomicnum==8:
                           count+=1
                    if count==2:
                        opbendvalue=round(1.5*71.94,2)
                    else:
                        opbendvalue=round(.05*71.94,2)
                else:
                    opbendvalue=round(poltype.defopbendval*71.94,2)

                defaultvalues.append(opbendvalue)
                firstprmline='opbend'+' '+str(poltypeclasses[0])+' '+str(poltypeclasses[1])+' '+'0'+' '+'0'+' '+str(opbendvalue)+'\n'
                newopbendprms.append(firstprmline)

    return newopbendprms,defaultvalues


def WriteOutList(poltype,ls,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    with open(filename, 'w') as filehandle:
        for listitem in ls:
            filehandle.write('%s\n' % listitem)

def ReadTorsionList(poltype,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newls=[]
    if os.stat(filename).st_size != 0:
        with open(filename, 'r') as filehandle:
            for line in filehandle:
                current = line[:-1]
                current=current.replace('[','').replace(']','')
                current=current.split(',')
                current=[int(i) for i in current]
                newls.append(current)

    return newls


def ReadTorTorList(poltype,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newls=[]
    if os.stat(filename).st_size != 0:
        with open(filename, 'r') as filehandle:
            for line in filehandle:
                current = line[:-1]
                current=current.replace('[','').replace(']','')
                current=current.split(',')
                current=[int(i) for i in current]
                newls.append(current)

    return newls



def ReadVdwList(poltype,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newls=[]
    if os.stat(filename).st_size != 0:
        with open(filename, 'r') as filehandle:
            for line in filehandle:
                current = line[:-1]
                newls.append(int(current))

    return newls


def CheckIfParametersExist(poltype,potentialmissingindices,prms):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    missingprmindices=[]
    for indices in potentialmissingindices:
        babelindices=[i+1 for i in indices]
        symtypes=[poltype.idxtosymclass[i] for i in babelindices]
        found=False
        symtypes=[str(i) for i in symtypes]
        for prmline in prms:
            linesplit=prmline.split()
            parms=[linesplit[1],linesplit[2]]
            if parms[0]==symtypes[0] and parms[1]==symtypes[1]:
               found=True
        if found==False:
            missingprmindices.append(indices)
    return missingprmindices
                
def WriteDictionaryToFile(poltype,dic,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    json.dump(dic, open(filename,'w'))        

def ReadDictionaryFromFile(poltype,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    return json.load(open(filename))

def GenerateSMARTSMatchLine(poltype,rdkitmol,rdkitindex):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    smarts=rdmolfiles.MolToSmarts(rdkitmol)
    smartsmol=Chem.MolFromSmarts(smarts)
    matches=rdkitmol.GetSubstructMatches(smartsmol)
    for match in matches:
        lastmatch=match
    lastindex=lastmatch.index(rdkitindex)
    atomorder=lastindex+1
    string='%'+' '+smarts+' % '+str(atomorder)
    return string

def GenerateAtomSMARTSMap(poltype,rdkitmol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    lines=[]
    descrip='"%s"'% poltype.molecprefix
    for atom in rdkitmol.GetAtoms():
        atomidx=atom.GetIdx()
        symbol=atom.GetSymbol()
        string=GenerateSMARTSMatchLine(poltype,rdkitmol,atomidx)
        line=' '+symbol+' '+descrip+' '+string+'\n'
        lines.append(line)
    return lines        

def AppendToSMARTSMapFile(poltype,lines,filename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(filename,'a')
    for line in lines:
        temp.write(line)
    temp.close()






def ReadExternalDatabase(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(poltype.externalparameterdatabase,'r')
    results=temp.readlines()
    temp.close()
    bondsmartsatomordertoparameters={} 
    anglesmartsatomordertoparameters={}
    strbndsmartsatomordertoparameters={}
    torsionsmartsatomordertoparameters={}
    tortorsmartsatomordertoparameters={}
    tortorsmartsatomordertogrid={}
    smartsatomordertotorvdwdb={}
    opbendsmartsatomordertoparameters={}
    vdwsmartsatomordertoparameters={}
    founddelim=False
    for line in results:
        linesplit=line.split()
        if len(linesplit)==0:
            continue
        if 'DELIM' in line:
            founddelim=True
        if linesplit[0]=='#':
            continue
        keyword=linesplit[0]
        linesplit=linesplit[1:]
        newline=' '.join(linesplit)
        linesplit=newline.split('%')
        smarts=linesplit[1].lstrip().rstrip()
        atomorderstring=linesplit[2]
        atomorderstringlist=atomorderstring.split(',')
        atomorderlist=tuple([int(i) for i in atomorderstringlist])
        prmstring=linesplit[3]
        prmstringlist=prmstring.split(',')
        prmlist=[float(i) for i in prmstringlist] 
        smartsatomorder=tuple([smarts,atomorderlist])
        if keyword=='bond':
            bondsmartsatomordertoparameters[smartsatomorder]=prmlist
        elif keyword=='angle':
            anglesmartsatomordertoparameters[smartsatomorder]=prmlist
        elif keyword=='strbnd':
            strbndsmartsatomordertoparameters[smartsatomorder]=prmlist
        elif keyword=='torsion':
            torsionsmartsatomordertoparameters[smartsatomorder]=prmlist
        elif keyword=='opbend':
            opbendsmartsatomordertoparameters[smartsatomorder]=prmlist
        elif keyword=='vdw':
            vdwsmartsatomordertoparameters[smartsatomorder]=prmlist
        elif keyword=='tortors':
            gridstring=linesplit[3].lstrip().rstrip()
            grid=gridstring.split()
            prmstring=linesplit[4]
            prmstringlist=prmstring.split(',')
            prmstringsplits=[s.lstrip().rstrip().split() for s in prmstringlist]
            tortorsmartsatomordertoparameters[smartsatomorder]=prmstringsplits
            tortorsmartsatomordertogrid[smartsatomorder]=grid
        if founddelim==True:
            smartsatomordertotorvdwdb[smartsatomorder]=True
        else:
            smartsatomordertotorvdwdb[smartsatomorder]=False
    return bondsmartsatomordertoparameters,anglesmartsatomordertoparameters,strbndsmartsatomordertoparameters,torsionsmartsatomordertoparameters,opbendsmartsatomordertoparameters,vdwsmartsatomordertoparameters,tortorsmartsatomordertoparameters,tortorsmartsatomordertogrid,smartsatomordertotorvdwdb


def ConvertToPoltypeClasses(poltype,torsionsmissing):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newtorsionsmissing=[]
    for sublist in torsionsmissing:
        newsublist=[i+1 for i in sublist]
        a,b,c,d=newsublist[:]
        batom=poltype.mol.GetAtom(b)
        catom=poltype.mol.GetAtom(c)
        aatom,datom = torgen.find_tor_restraint_idx(poltype,poltype.mol,batom,catom)
        newa=aatom.GetIdx()
        newd=datom.GetIdx()
        sorttor=torfit.sorttorsion(poltype,[poltype.idxtosymclass[a],poltype.idxtosymclass[b],poltype.idxtosymclass[c],poltype.idxtosymclass[d]])
        newsorttor=torfit.sorttorsion(poltype,[poltype.idxtosymclass[newa],poltype.idxtosymclass[b],poltype.idxtosymclass[c],poltype.idxtosymclass[newd]])

        if sorttor not in newtorsionsmissing:
            newtorsionsmissing.append(sorttor)
        if newsorttor not in newtorsionsmissing:
            newtorsionsmissing.append(newsorttor)

    return newtorsionsmissing



def GrabKeysFromValue(poltype,dic,thevalue):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    keylist=[]
    for key,value in dic.items():
        if value==thevalue:
            keylist.append(key)
    return keylist



def GrabAllPossibleMoleculeIndices(poltype,typeindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    thelist=[]
    for typenum in typeindices:
        allindices=GrabKeysFromValue(poltype,poltype.idxtosymclass,typenum)
        thelist.append(allindices)
    listofmoleculeindices=list(itertools.product(*thelist))


    return listofmoleculeindices


def MatchExternalSMARTSToMolecule(poltype,rdkitmol,smartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb,torsionindicestoextsmarts=None):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indicestoextsmartsmatchlength={}
    indicestoextsmarts={}
    indicestoextsmartsatomorder={}
    indicestosmartslist={}
    indicestomatchlist={}
    indicestoatomorderlist={}
    indicestosmartsatomorder={}
    indicestoprmlist={}
    restrictedsmarts=[]
    failedmatches=[]
    if poltype.quickdatabasesearch==False:
        if torsionindicestoextsmarts!=None:
            for torsionindices,extsmarts in torsionindicestoextsmarts.items():
                if extsmarts not in restrictedsmarts:
                    restrictedsmarts.append(extsmarts)
        listsmartsatomordertoparameters=list(smartsatomordertoparameters.keys())
        for k in tqdm(range(len(listsmartsatomordertoparameters)),desc='Searching SMARTS database'):
            smartsatomorder=listsmartsatomordertoparameters[k]
            parameters=smartsatomordertoparameters[smartsatomorder]
            smarts=smartsatomorder[0]
            fromtorvdwdb=smartsatomordertotorvdwdb[smartsatomorder]
            atomorderlist=smartsatomorder[1]
            if len(restrictedsmarts)!=0:
                if smarts not in restrictedsmarts and fromtorvdwdb==True:
                    continue
            if (len(restrictedsmarts)==0 and len(atomorderlist)==1)and fromtorvdwdb==True:
                continue 
            substructure = Chem.MolFromSmarts(smarts)
            mols = [rdkitmol,substructure]
            res=rdFMCS.FindMCS(mols)
            atomnum=res.numAtoms
            smartsmcs=res.smartsString
            diditmatch=False
            diditmatchexactly=rdkitmol.HasSubstructMatch(substructure)
            if fromtorvdwdb==False and diditmatchexactly==False:
                continue
            if atomnum>=len(atomorderlist):

                mcssubstructure = Chem.MolFromSmarts(smartsmcs)
                mcsmatches=substructure.GetSubstructMatches(mcssubstructure)
                if len(mcsmatches)>0:
                    mcsmatch=mcsmatches[0]
                    mcsindices=list(range(len(mcsmatch)))
                    mcssmartsindextosmartsmolindex=dict(zip(mcsindices,mcsmatch)) 
                    smartsmolindexlist=[i-1 for i in atomorderlist]
                    foundindices=True
                    
                    for i in smartsmolindexlist:
                        if i not in mcssmartsindextosmartsmolindex.values():
                            foundindices=False
                    if foundindices==True:
                        matches=rdkitmol.GetSubstructMatches(mcssubstructure)
                        if len(matches)>0:
                            thematch=matches[0]
                            theindices=list(range(len(thematch)))
                            mcssmartsindextomolindex=dict(zip(theindices,thematch)) 
                            smartsmolindextomcssmartsindex={v: k for k, v in mcssmartsindextosmartsmolindex.items()} 
                            mcssmartsindices=[smartsmolindextomcssmartsindex[i] for i in smartsmolindexlist]
                            indices=[mcssmartsindextomolindex[i] for i in mcssmartsindices]
                            mcsorder=[i+1 for i in indices]
                            
                            for match in matches:
                                allin=True
                                for idx in indices:
                                    if idx not in match:
                                        allin=False
                                if allin==True:
                                    matcharray=match
                                    break
                            nindexes=[]
                            for idx in indices: 
                                neighborindexes=indextoneighbidxs[idx]
                                for nidx in neighborindexes:
                                    if nidx not in nindexes:
                                        nindexes.append(nidx)
                                    if len(indices)==1 and len(neighborindexes)==1:
                                        nneighborindexes=indextoneighbidxs[nidx]
                                        for nnidx in nneighborindexes:
                                            if nnidx not in nindexes:
                                                nindexes.append(nnidx)
            
                            diditmatch=CheckIfNeighborsExistInSMARTMatch(poltype,nindexes,matcharray)  
                            
                            if fromtorvdwdb==False:
                                diditmatch=True # only care about neighbor matching for large SMARTS from part of DB genreated by fragmenter, top of DB can have shorter SMARTS 
                            if diditmatch==True:
                                moleculeindices=tuple(indices)
                                if len(indices)==4: # then torsion matches need to be consistent for all torsion around bond
                                    b=indices[1]
                                    c=indices[2]
                                    batombabel=poltype.mol.GetAtom(b+1)
                                    catombabel=poltype.mol.GetAtom(c+1)
                                    bisinring=batombabel.IsInRing()
                                    cisinring=catombabel.IsInRing()
                                    bhyb=batombabel.GetHyb()
                                    chyb=catombabel.GetHyb()
                                    if (bisinring==True and cisinring==True) and (bhyb==2 and chyb==2):
                                        continue
                                if moleculeindices not in indicestosmartslist.keys():
                                    indicestosmartslist[moleculeindices]=[]
                                    indicestomatchlist[moleculeindices]=[]
                                    indicestoatomorderlist[moleculeindices]=[]
                                    indicestoprmlist[moleculeindices]=[]
                                    indicestosmartsatomorder[moleculeindices]=[]

                                if smartsmcs not in indicestosmartslist[moleculeindices]:
                                    mcssmartsatomorder=tuple([smartsmcs,tuple(mcsorder)])
                                    indicestosmartslist[moleculeindices].append(smartsmcs)
                                    indicestomatchlist[moleculeindices].append(matcharray)
                                    indicestoatomorderlist[moleculeindices].append(mcssmartsatomorder)
                                    indicestosmartsatomorder[moleculeindices].append(smartsatomorder)
                                    indicestoprmlist[moleculeindices].append(parameters)
            else:
                failedmatches.append(smarts)
        rotatablebondtosmartslist={}
        for moleculeindices,smartslist in indicestosmartslist.items():
            if len(moleculeindices)==4: 
                b=moleculeindices[1]
                c=moleculeindices[2]
                if b<c:
                    rotbnd=[b,c]
                elif b>c:
                    rotbnd=[c,b]
                else:
                    rotbnd=[b,c] 
                rotbnd=tuple(rotbnd)
                if rotbnd not in rotatablebondtosmartslist.keys():
                    rotatablebondtosmartslist[rotbnd]=[]
                for smarts in smartslist:
                    if smarts not in rotatablebondtosmartslist[rotbnd]:
                        rotatablebondtosmartslist[rotbnd].append(smarts)
        rotatablebondtotruesmartslist={}
        for rotbnd,allsmartslist in rotatablebondtosmartslist.items():
            for smarts in allsmartslist:
                allin=True
                for moleculeindices,smartslist in indicestosmartslist.items():
                    if len(moleculeindices)==4: 
                       b=moleculeindices[1]
                       c=moleculeindices[2]
                       if b<c:
                           otherrotbnd=[b,c]
                       elif b>c:
                           otherrotbnd=[c,b]
                       else:
                           otherrotbnd=[b,c] 
                       otherrotbnd=tuple(otherrotbnd)
                       if rotbnd==otherrotbnd: 
                           if smarts not in smartslist:
                               allin=False
                if allin==True: 
                    if rotbnd not in rotatablebondtotruesmartslist.keys():
                        rotatablebondtotruesmartslist[rotbnd]=[]
                    if smarts not in rotatablebondtotruesmartslist[rotbnd]:
                        rotatablebondtotruesmartslist[rotbnd].append(smarts)

        rotatablebondtotruesmarts={}
        for rotbnd,smartslist in rotatablebondtotruesmartslist.items():
            smartslistlen=[len(i) for i in smartslist]
            maxlen=max(smartslistlen)
            maxidx=smartslistlen.index(maxlen)
            maxsmarts=smartslist[maxidx]
            rotatablebondtotruesmarts[rotbnd]=maxsmarts

        for moleculeindices,smartslist in indicestosmartslist.items():
           if len(moleculeindices)==4: 
               b=moleculeindices[1]
               c=moleculeindices[2]
               if b<c:
                   rotbnd=[b,c]
               elif b>c:
                   rotbnd=[c,b]
               else:
                   rotbnd=[b,c] 
               rotbnd=tuple(rotbnd)
               if rotbnd in rotatablebondtotruesmarts.keys():
                   maxsmarts=rotatablebondtotruesmarts[rotbnd]
               else:
                   smartslistlen=[len(i) for i in smartslist]
                   maxlen=max(smartslistlen)
                   maxidx=smartslistlen.index(maxlen)
                   maxsmarts=smartslist[maxidx]
          
        
           else:
               smartslistlen=[len(i) for i in smartslist]
               maxlen=max(smartslistlen)
               maxidx=smartslistlen.index(maxlen)
               maxsmarts=smartslist[maxidx]
           for j in range(len(indicestosmartsatomorder[moleculeindices])):
               originalsmartsatomorder=indicestosmartsatomorder[moleculeindices][j]
               thesmarts=indicestosmartslist[moleculeindices][j]
               if thesmarts==maxsmarts: 
                   break


           smarts=maxsmarts
           matcharray=indicestomatchlist[moleculeindices][j]
           smartsatomorder=indicestoatomorderlist[moleculeindices][j]
           prms=indicestoprmlist[moleculeindices][j]
           babelindices=[i+1 for i in moleculeindices]
           typeindices=[poltype.idxtosymclass[i] for i in babelindices]
           listofmoleculeindices=GrabAllPossibleMoleculeIndices(poltype,typeindices) 
           for babelindices in listofmoleculeindices:
               moleculeindices=tuple([i-1 for i in babelindices])
               if moleculeindices not in indicestoextsmarts.keys() and moleculeindices[::-1] not in indicestoextsmarts.keys():
                   indicestoextsmartsmatchlength[moleculeindices]=len(matcharray)
                   indicestoextsmarts[moleculeindices]=smarts
                   indicestoextsmartsatomorder[moleculeindices]=smartsatomorder
                   smartsatomordertoparameters[smartsatomorder]=prms
    return indicestoextsmartsmatchlength,indicestoextsmarts,indicestoextsmartsatomorder,smartsatomordertoparameters



def CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,indicestoextsmartsmatchlength,indicesforprmtoparametersmarts,indicesforprmtosmarts,indicestoextsmarts,indicesforprmtomatchallneighbs,indicestoextsmartsatomorders):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newindicestoextsmarts={}
    newindicestoextsmartsatomorder={}
    for indices,extsmartsmatchlength in indicestoextsmartsmatchlength.items():
        if indices in indicesforprmtoparametersmarts.keys():
            smartsls=indicesforprmtosmarts[indices]
            smarts=smartsls[0]
            matchallneighbs=indicesforprmtomatchallneighbs[indices]
            if indices in indicestoextsmartsatomorders.keys():
                smartsatomorder=indicestoextsmartsatomorders[indices]
                extsmarts=indicestoextsmarts[indices]
            else:
                smartsatomorder=indicestoextsmartsatomorders[indices[::-1]]
                extsmarts=indicestoextsmarts[indices[::-1]]
            poormatch=True
            substructure = Chem.MolFromSmarts(smarts)
            substructurenumatoms=substructure.GetNumAtoms()
            if extsmartsmatchlength>substructurenumatoms or poormatch==True:
                del indicesforprmtosmarts[indices]
                del indicesforprmtoparametersmarts[indices]
                newindicestoextsmarts[indices]=extsmarts
                newindicestoextsmartsatomorder[indices]=smartsatomorder
        elif indices[::-1] in indicesforprmtoparametersmarts.keys():
            smartsls=indicesforprmtosmarts[indices[::-1]]
            smarts=smartsls[0]
            matchallneighbs=indicesforprmtomatchallneighbs[indices[::-1]]
            if indices[::-1] in indicestoextsmartsatomorders.keys():
                smartsatomorder=indicestoextsmartsatomorders[indices[::-1]]
                extsmarts=indicestoextsmarts[indices[::-1]]
            else:
                smartsatomorder=indicestoextsmartsatomorders[indices]
                extsmarts=indicestoextsmarts[indices]

            poormatch=True
            substructure = Chem.MolFromSmarts(smarts)
            substructurenumatoms=substructure.GetNumAtoms()
            if extsmartsmatchlength>substructurenumatoms:
                del indicesforprmtosmarts[indices[::-1]]
                del indicesforprmtoparametersmarts[indices[::-1]]
                 
                newindicestoextsmarts[indices[::-1]]=extsmarts
                newindicestoextsmartsatomorder[indices[::-1]]=smartsatomorder
    return indicesforprmtoparametersmarts,indicesforprmtosmarts,newindicestoextsmarts,newindicestoextsmartsatomorder

def AddExternalDatabaseSMARTSMatchParameters(poltype,prms,indicestoextsmarts,smartsatomordertoparameters,keyword,indicestoextsmartsatomorder,smartsatomordertogrid=None):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    poltypeclassestosmartsatomordersext={}
    for indices,extsmarts in indicestoextsmarts.items():
        thesmartsatomorder=indicestoextsmartsatomorder[indices]
        theatomorder=thesmartsatomorder[1]
        for smartsatomorder,parameters in smartsatomordertoparameters.items():
            smarts=smartsatomorder[0]
            atomorder=smartsatomorder[1]
            
            if smarts==extsmarts and atomorder==theatomorder:
                line=keyword+' '
                babelindices=[i+1 for i in indices]
                poltypeclasses=[poltype.idxtosymclass[i] for i in babelindices]
                for poltypeclass in poltypeclasses:
                    line+=str(poltypeclass)+' '
                if keyword=='opbend':
                    line+='0 0 '
                if keyword=='tortors':
                    grid=smartsatomordertogrid[smartsatomorder]
                    gridstring=' '.join(grid)
                    line+=gridstring
                    line+='\n'
                    prms.append(line)
                    for point in parameters:
                        string=' '.join(point)
                        string+='\n'
                        prms.append(string)
                    pass
                elif keyword=='torsion':
                    for prmidx in range(len(parameters)):
                        prm=parameters[prmidx]
                        if prmidx==1:
                            phase=180
                        else:
                            phase=0
                        fold=prmidx+1
                        line+=str(prm)+' '+str(phase)+' '+str(fold)+' ' 
                    line+='\n' 

                else:
               
                    for prm in parameters:
                        line+=str(prm)+' '
                    line+='\n' 
                prms.append(line)

                poltypeclassestosmartsatomordersext[tuple([tuple(poltypeclasses)])]=smartsatomorder
    return prms,poltypeclassestosmartsatomordersext


def FilterBondSMARTSEnviorment(poltype,bondindicestosmartsatomorders,bondindicestotinkerclasses):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    opbendbondindicestosmartsatomorders={}
    opbendbondindicestotinkerclasses={}
    for bondindices,smartsatomorders in bondindicestosmartsatomorders.items(): 
        smarts=smartsatomorders[0]
        brackcount=CountBrackets(poltype,smarts)
        revbondindices=bondindices[::-1]
        if brackcount==2:
            pass
        else:
            temp=bondindicestotinkerclasses[bondindices]
            revtemp=temp[::-1]
            opbendbondindicestotinkerclasses[bondindices]=temp
            opbendbondindicestotinkerclasses[revbondindices]=revtemp
            opbendbondindicestosmartsatomorders[bondindices]=smartsatomorders 
            opbendbondindicestosmartsatomorders[revbondindices]=smartsatomorders 

    return opbendbondindicestotinkerclasses,opbendbondindicestosmartsatomorders

def CountBrackets(poltype,string):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    count=0
    for e in string:
        if e=='[':
            count+=1
    return count

def CheckTrigonalCenters(poltype,listofbondsforprm,mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    opbendbondindicestotrigonalcenterbools={}
    for bondindices in listofbondsforprm:
        bondindices=tuple(bondindices)
        boolarray=[]
        babelindices=[i+1 for i in bondindices]
        atoms=[mol.GetAtom(i) for i in babelindices]     
        for a in atoms:
            hyb=a.GetHyb()
            numhyds=0
            neighbs=list(openbabel.OBAtomAtomIter(a))
            for natom in neighbs:
                natomicnum=natom.GetAtomicNum()
                if natomicnum==1:
                    numhyds+=1
            idx=a.GetIdx()
            atomicnum=a.GetAtomicNum()
            if len(neighbs)==3 and hyb==2:
                if atomicnum==7 and numhyds==2:
                    boolarray.append(False)

                else:
                    boolarray.append(True)
            else:
                boolarray.append(False)
        opbendbondindicestotrigonalcenterbools[bondindices]=boolarray
        revboolarray=boolarray[::-1]
        revbondindices=bondindices[::-1]
        opbendbondindicestotrigonalcenterbools[revbondindices]=revboolarray
    return opbendbondindicestotrigonalcenterbools


def CorrectPitorEnergy(poltype,torsionprms,torsiontopitor):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newtorsionprms=[]
    torsiontocount={}
    middletocount={}
     
    for torsion in torsiontopitor.keys():
        middle=tuple([torsion[1],torsion[2]])
        if middle not in middletocount.keys():
            middletocount[middle]=0
        middletocount[middle]+=1
    for torsion in torsiontopitor.keys():
        middle=tuple([torsion[1],torsion[2]])
        count=middletocount[middle]
        torsiontocount[torsion]=count
    for torline in torsionprms:
        torlinesplit=torline.split()
        tor=tuple([int(torlinesplit[1]),int(torlinesplit[2]),int(torlinesplit[3]),int(torlinesplit[4])])
        if tor in torsiontopitor:
            if tor in torsiontopitor.keys():
                pitorline=torsiontopitor[tor]
            else:
                pitorline=torsiontopitor[tor[::-1]]

            count=torsiontocount[tor]
            pitorlinesplit=pitorline.split()
            prm=float(pitorlinesplit[3])/count
            torprm=float(torlinesplit[8])
            newtorprm=prm+torprm
            torlinesplit[8]=str(newtorprm)
        torline=' '.join(torlinesplit)+'\n'
        newtorsionprms.append(torline)
    return newtorsionprms


def FindMissingParameters(poltype,indicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    missing=[]
    for indices,smartsatomorders in indicestosmartsatomorders.items():
        nindexes=[]
        for idx in indices: 
            neighborindexes=indextoneighbidxs[idx]
            for nidx in neighborindexes:
                if nidx not in nindexes:
                    nindexes.append(nidx)
                if len(indices)==1 and len(neighborindexes)==1:
                    nneighborindexes=indextoneighbidxs[nidx]
                    for nnidx in nneighborindexes:
                        if nnidx not in nindexes:
                            nindexes.append(nnidx) 
        smarts=smartsatomorders[0]
        substructure = Chem.MolFromSmarts(smarts)
        matches=rdkitmol.GetSubstructMatches(substructure)
        for match in matches:
            allin=True
            for idx in indices:
                if idx not in match:
                    allin=False
            if allin==True:
                matcharray=match
                break

         
        check=CheckIfNeighborsExistInSMARTMatch(poltype,nindexes,matcharray)
        if check==False or '*' in smarts or '~' in smarts:
            if len(indices)==1: # vdw
                index=indices[0]
                if poltype.onlyvdwatomindex!=None:
                    idx=indices[0]+1
                    if idx==poltype.onlyvdwatomindex:
                        missing.append(indices)
                else:
                    missing.append(indices)
            else:
                missing.append(indices)
        else:
            if len(indices)==1: # vdw
                index=indices[0]
                if poltype.onlyvdwatomindex==index:
                    missing.append(indices)





    return missing

def ReduceMissingVdwByTypes(poltype,vdwmissing):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    reducedvdwmissing=[]
    typesfound=[]
    for vdwarray in vdwmissing:
        vdwatomidx=vdwarray[0]
        vdwbabelidx=vdwatomidx+1
        vdwtype=poltype.idxtosymclass[vdwbabelidx]
        if vdwtype in typesfound:
            continue
        else:
            typesfound.append(vdwtype)
            reducedvdwmissing.append(vdwbabelidx)
    return reducedvdwmissing



def ConvertToBabelList(poltype,listforprm):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    babellist=[]
    for ls in listforprm:
        babells=[i+1 for i in ls]
        babellist.append(babells)

    return babellist



def AddReverseKeys(poltype,tinkerclassestopoltypeclasses):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newtinkerclassestopoltypeclasses={}
    for tinkerclasses,poltypeclasses in tinkerclassestopoltypeclasses.items():
        revtinkerclasses=tinkerclasses[::-1]
        revpoltypeclasses=[]
        for ls in poltypeclasses:
            revls=ls[::-1]
            revpoltypeclasses.append(revls) 
        newtinkerclassestopoltypeclasses[tinkerclasses]=poltypeclasses
        newtinkerclassestopoltypeclasses[revtinkerclasses]=revpoltypeclasses
    return newtinkerclassestopoltypeclasses


def TinkerClassesToTrigonalCenter(poltype,opbendbondindicestotinkerclasses,opbendbondindicestotrigonalcenterbools):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    opbendtinkerclassestotrigonalcenterbools={}
    for bondindices,tinkerclasses in opbendbondindicestotinkerclasses.items():
        boolarray=opbendbondindicestotrigonalcenterbools[bondindices]
        opbendtinkerclassestotrigonalcenterbools[tuple(tinkerclasses)]=boolarray
        revtinkerclasses=tinkerclasses[::-1]
        revboolarray=boolarray[::-1]
        opbendtinkerclassestotrigonalcenterbools[tuple(revtinkerclasses)]=revboolarray

    return opbendtinkerclassestotrigonalcenterbools


def FilterDictionaries(poltype,dics,ls):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newdics=[]
    for dic in dics:
        newdic={}
        for key,value in dic.items():
            if key in ls:
                newdic[key]=value
        newdics.append(newdic)
    return newdics

def ConvertListOfListToListOfTuples(poltype,listoflist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    listoftuples=[]
    for item in listoflist:
        tup=tuple(item)
        listoftuples.append(tup)
    return listoftuples


def AddExternalDatabaseMatches(poltype, indicestosmartsatomorder,extindicestoextsmarts,smartsatomordertoparameters):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newindicestosmartsatomorder=indicestosmartsatomorder.copy()
    for smartsatomorder in smartsatomordertoparameters.keys():
        smarts=smartsatomorder[0]
        for indices,extsmarts in extindicestoextsmarts.items():
            if smarts==extsmarts:
                newindicestosmartsatomorder[indices]=list(smartsatomorder)         
        
    return newindicestosmartsatomorder   
   
def AngleGuess(poltype,ita,itb,itc,iva,ivb,ivc):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    radian=57.29577951
    angunit = 1.0 / radian**2
    if (itb==6):
       if (ita==1):
          if (ivb==4):
             if (itc==1):
                angguess = 34.50
             elif (itc==6):
                angguess = 38.0
             elif (itc==7):
                angguess = 50.60
             elif (itc==8):
                angguess = 51.50
             elif (itc==9):
                angguess = 50.0
             else:
                angguess = 35.0
          elif (ivb==3):
             angguess = 32.00
          else:
             angguess = 32.00
       elif (ita==6):
          if (ivb==4):
             if (itc==6):
                angguess = 60.00
             elif (itc==7):
                angguess = 80.00
             elif (itc==8):
                angguess = 88.00
             elif (itc==9):
                angguess = 89.00
             elif (itc==14):
                angguess = 65.00
             elif (itc==15):
                angguess = 60.00
             elif (itc==16):
                angguess = 53.20
             elif(itc==17):
                angguess = 55.00
             else:
                angguess = 50.00
          elif (ivb==3):
             angguess = 60.00
          else:
             angguess = 60.00
       elif (ita==8):
          if (ivb==4):
             if (itc==8):
                angguess = 65.00
             elif (itc==9):
                angguess = 65.00
             elif (itc==15):
                angguess = 60.00
             elif (itc==16):
                angguess = 65.00
             else:
                angguess = 65.00
             
          elif (ivb==3):
             angguess = 50.00
          else:
             angguess = 60.00
          
       else:
          angguess = 60.00
       
    elif (itb==8):
       if (ita==1):
          if (itc==1):
             angguess = 34.05
          elif (itc==6):
             angguess = 65.00
          else:
             angguess = 60.00
          
       elif (ita==6):
          if (itc==6):
             angguess = 88.50
          elif (itc==8):
             if (iva==1 or ivc==1):
                angguess = 122.30
             else:
                angguess = 85.00
             
          elif (itc==15):
             angguess = 80.30
          else:
             angguess = 80.0
          
       else:
          angguess = 80.0
       
    elif (itb==15):
       if (ita==1):
          angguess = 30.0
       elif (ita==6):
          if (itc==6):
             angguess = 75.00
          elif (itc==8):
             angguess = 80.00
          else:
             angguess = 75.00
         
       elif (ita==8):
          if (itc==8):
             if (iva==1 and ivc==1):
                angguess = 89.88
             elif (iva==1 or ivc==1):
                angguess = 75.86
             else:
                angguess = 65.58
             
          else:
             angguess = 70.00
          
       else:
          angguess = 75.00
       
    elif (itb==16):
       if (ita==1):
          angguess = 30.00
       elif (ita==6):
          if (itc==16):
             angguess = 72.00
          else:
             angguess = 80.00
          
       elif (ita==8):
          if (itc==8):
             if (iva==1 and ivc==1):
                angguess = 168.00
             elif (iva==1 or ivc==1):
                angguess = 85.00
             else:
                angguess = 80.00
             
          elif (itc==16):
             angguess = 75.00
          else:
             angguess = 75.00
          
       else:
          angguess = 75.00
       
    elif (ita==1):
       angguess = 35.00
    else:
       angguess = 65.00
    
    angguess = angguess / (angunit*radian**2)
    return angguess
 
 
def BondGuess(poltype,ita,itb,iva,ivb):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    bndunit=1
    if (ita==1):
         if (itb==6):
            if (ivb==3):
               bndguess = 410.0
            elif (ivb==4):
               bndguess = 400.0
            else:
               bndguess = 400.0
            
         elif (itb==7):
            bndguess = 520.0
         elif (itb==8):
            bndguess = 560.0
         elif (itb==9):
            bndguess = 500.0
         elif (itb==14):
            bndguess = 200.0
         elif (itb==15):
            bndguess = 230.0
         elif (itb==16):
            bndguess = 260.0
         else:
            bndguess = 300.0
         
    elif (ita==6):
       if (itb==6):
          if (iva==3 and ivb==3):
             bndguess = 680.0
          elif (iva==4 or ivb==4):
             bndguess = 385.0
          else:
             bndguess = 350.0
          
       elif (itb==7):
          if (iva==3 and ivb==2):
             bndguess = 435.0
          elif (iva==3 and ivb==3):
             bndguess = 250.0
          elif (iva==4):
             bndguess = 400.0
          else:
             bndguess = 450.0
          
       elif (itb==8):
          if (ivb==1):
             bndguess = 680.0
          elif (ivb==2):
             bndguess = 465.0
          else:
             bndguess = 465.0
          
       elif (itb==9):
          bndguess = 350.0
       elif (itb==14):
          bndguess = 350.0
       elif (itb==15):
          bndguess = 350.0
       elif (itb==16):
          bndguess = 216.0
       elif (itb==17):
          bndguess = 350.0
       else:
          bndguess = 450.0
       
    elif (ita==7):
       if (itb==7):
          if (iva==1):
             bndguess = 1613.0
          elif (iva==2 and ivb==2):
             bndguess = 950.0
          else:
             bndguess = 850.0
          
       elif (itb==8):
          if (ivb==1):
             bndguess = 900.0
          else:
             bndguess = 750.0
          
       elif (itb==14):
          bndguess = 450.0
       elif (itb==15):
          bndguess = 500.0
       elif (itb==16):
          bndguess = 550.0
       else:
          bndguess = 600.0
       
    elif (ita==8):
       if (itb==8):
          bndguess = 750.0
       elif (itb==14):
          bndguess = 500.0
       elif (itb==15):
          if (iva==2):
             bndguess = 450.0
          elif (iva==1):
             bndguess = 775.0
          else:
             bndguess = 450.0
          
       elif (itb==16):
          bndguess = 606.0
       elif (itb==17):
          bndguess = 500.0
       else:
          bndguess = 600.0
       
    elif (ita==14):
       if (itb==14):
          bndguess = 400.0
       elif (itb==15):
          bndguess = 450.0
       elif (itb==16):
          bndguess = 500.0
       elif (itb==17):
          bndguess = 650.0
       else:
          bndguess = 450.0
      
    elif (ita==16):
       if (itb==16):
          bndguess = 188.0
       else:
          bndguess = 250.0
       
    elif (ita==17):
       bndguess = 300.0
    else:
       bndguess = 350.0
      
    bndguess = bndguess / bndunit
    return bndguess 

def ReverseDictionaryValueList(poltype,keytovalues):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    valuetokey={}
    for key,value in keytovalues.items():
        for ls in value:
            valuetokey[tuple(ls)]=list(key)
    return valuetokey

def ReverseDictionary(poltype,keytovalues):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    valuetokey={}
    for key,value in keytovalues.items():
        valuetokey[tuple(value)]=list(key)
    return valuetokey

def RepresentsInt(s):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    try: 
        int(s)
        return True
    except ValueError:
        return False

def ReadDatabaseSmartsMap(poltype,databasepath): # smartsString classNumber className # comments
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    smartstoatomclass={}
    atomclasstoclassname={}
    atomclasstocomment={}
    temp=open(databasepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)>0:
            first=linesplit[0]
            if '#' not in first:
                smarts=linesplit[0]
                if RepresentsInt(linesplit[1])==True:
                    tinkerclass=int(linesplit[1])
                else:
                    tinkerclass=linesplit[1]
                tinkername=linesplit[2]
                comment=' '.join(linesplit[3:])
                comment=comment.replace('\n','').replace('#','').lstrip().rstrip()
                smartstoatomclass[smarts]=tinkerclass
                atomclasstocomment[tinkerclass]=comment
                atomclasstoclassname[tinkerclass]=tinkername
             

    return smartstoatomclass, atomclasstoclassname, atomclasstocomment
 
def ReadDatabaseSmartsMapPolarize(poltype,databasepath): # smartsString className # comments
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    smartstoatomclass={}
    atomclasstocomment={}
    smartstoactualcomment={}
    temp=open(databasepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)>0:
            first=linesplit[0]
            if '# ' in first or first=='#':
                continue
            else:
                smarts=linesplit[0]
                if RepresentsInt(linesplit[1])==True:
                    tinkerclass=int(linesplit[1])
                else:
                    tinkerclass=linesplit[1]
                comment=' '.join(linesplit[3:])
                comment=comment.replace('\n','').replace('#','').lstrip().rstrip()
                smartstoatomclass[smarts]=tinkerclass
                atomclasstocomment[tinkerclass]=comment
                smartstoactualcomment[smarts]=comment
             

    return smartstoatomclass, atomclasstocomment,smartstoactualcomment
 
def MatchAllSmartsToAtomIndices(poltype,smartstoatomclass): #rdkit 0 index based
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atomindextoallsmarts={}
    atomindextoallsmartsmatches={}
    for smarts in smartstoatomclass.keys():
        substructure = Chem.MolFromSmarts(smarts)
        atomclass=smartstoatomclass[smarts]
        diditmatch=poltype.rdkitmol.HasSubstructMatch(substructure)
        if diditmatch==True:
            matches=list(poltype.rdkitmol.GetSubstructMatches(substructure))
            newmatches=[]
            for match in matches:
                newmatches.append(match)
                lastatom=match[-1]
                firstatom=match[0]
                firsttype=poltype.idxtosymclass[firstatom+1]
                secondtype=poltype.idxtosymclass[lastatom+1]
                if firsttype==secondtype:
                    newmatches.append(match[::-1])
                otheratomindices=[firstatom]
                for idx,typenum in poltype.idxtosymclass.items():
                    rdkitidx=idx-1
                    if typenum==firsttype and rdkitidx!=firstatom:
                        otheratomindices.append(rdkitidx)

            for match in newmatches:
                for atomindex in otheratomindices:
                    if atomindex not in atomindextoallsmarts.keys():
                        atomindextoallsmarts[atomindex]=[]
                        atomindextoallsmartsmatches[atomindex]=[] 
                    if smarts not in atomindextoallsmarts[atomindex]:
                        atomindextoallsmarts[atomindex].append(smarts) 
                    if match not in atomindextoallsmartsmatches[atomindex]:
                        atomindextoallsmartsmatches[atomindex].append(match)   
               
    return atomindextoallsmarts,atomindextoallsmartsmatches


def MapIndicesToCommentsAtom(poltype,atomindextoallsmarts,smartstocomment,listofatomsforprm):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atomcommentstolistofsmartslist={}
    atomindicestolistofatomcomments={}
    for atoms in listofatomsforprm:
        aindex=atoms[0]
        if aindex in atomindextoallsmarts.keys():
            asmartslist=atomindextoallsmarts[aindex]
            combs = list(itertools.product(asmartslist)) 
            for comb in combs:
                asmarts=comb[0]
                acomment=smartstocomment[asmarts]
                comments=tuple([acomment])
                smartslist=[asmarts]
                if comments not in atomcommentstolistofsmartslist.keys():
                    atomcommentstolistofsmartslist[comments]=[]
                if tuple(atoms) not in atomindicestolistofatomcomments.keys(): 
                    atomindicestolistofatomcomments[tuple(atoms)]=[]
                if smartslist not in atomcommentstolistofsmartslist[comments]: 
                    atomcommentstolistofsmartslist[comments].append(smartslist)   
                if comments not in atomindicestolistofatomcomments[tuple(atoms)]: 
                    atomindicestolistofatomcomments[tuple(atoms)].append(comments)    
    return atomcommentstolistofsmartslist,atomindicestolistofatomcomments


def MapIndicesToCommentsBondAngle(poltype,atomindextoallsmarts,smartstocomment,listofbondsforprm,listofanglesforprm):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    bondcommentstolistofsmartslist={}
    bondindicestolistofbondcomments={}
    for bond in listofbondsforprm:
        aindex=bond[0] 
        bindex=bond[1] 
        asmartslist=atomindextoallsmarts[aindex]
        bsmartslist=atomindextoallsmarts[bindex]
        combs = list(itertools.product(asmartslist,bsmartslist)) 
        for comb in combs:
            asmarts=comb[0]
            bsmarts=comb[1]
            acomment=smartstocomment[asmarts]
            bcomment=smartstocomment[bsmarts]
            comments=tuple([acomment,bcomment])
            smartslist=[asmarts,bsmarts]
            if comments not in bondcommentstolistofsmartslist.keys():
                bondcommentstolistofsmartslist[comments]=[]
            if tuple(bond) not in bondindicestolistofbondcomments.keys(): 
                bondindicestolistofbondcomments[tuple(bond)]=[]
            if smartslist not in bondcommentstolistofsmartslist[comments]: 
                bondcommentstolistofsmartslist[comments].append(smartslist)   
            if comments not in bondindicestolistofbondcomments[tuple(bond)]: 
                bondindicestolistofbondcomments[tuple(bond)].append(comments)    

    anglecommentstolistofsmartslist={}
    angleindicestolistofanglecomments={}
    for angle in listofanglesforprm:
        aindex=angle[0] 
        bindex=angle[1] 
        cindex=angle[2] 
        asmartslist=atomindextoallsmarts[aindex]
        bsmartslist=atomindextoallsmarts[bindex]
        csmartslist=atomindextoallsmarts[cindex]
        combs = list(itertools.product(asmartslist,bsmartslist,csmartslist)) 
        for comb in combs:
            asmarts=comb[0]
            bsmarts=comb[1]
            csmarts=comb[2]
            acomment=smartstocomment[asmarts]
            bcomment=smartstocomment[bsmarts]
            ccomment=smartstocomment[csmarts]
            comments=tuple([acomment,bcomment,ccomment])
            smartslist=[asmarts,bsmarts,csmarts]
            if comments not in anglecommentstolistofsmartslist.keys():
                anglecommentstolistofsmartslist[comments]=[]
            if tuple(angle) not in angleindicestolistofanglecomments.keys(): 
                angleindicestolistofanglecomments[tuple(angle)]=[]
            if smartslist not in anglecommentstolistofsmartslist[comments]: 
                anglecommentstolistofsmartslist[comments].append(smartslist)   
            if comments not in angleindicestolistofanglecomments[tuple(angle)]: 
                angleindicestolistofanglecomments[tuple(angle)].append(comments)    

    return bondcommentstolistofsmartslist,bondindicestolistofbondcomments,anglecommentstolistofsmartslist,angleindicestolistofanglecomments




def MapIndicesToClasses(poltype,atomindextoallsmarts,smartstoatomclass,listofbondsforprm,listofanglesforprm,planarbonds):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    bondclassestolistofsmartslist={}
    angleclassestolistofsmartslist={}
    strbndclassestolistofsmartslist={}
    opbendclassestolistofsmartslist={}
    bondindicestolistofbondclasses={}
    angleindicestolistofangleclasses={}
    strbndindicestolistofstrbndclasses={}
    opbendindicestolistofopbendclasses={}
    for bond in listofbondsforprm:
        aindex=bond[0]
        bindex=bond[1] 
        if aindex in atomindextoallsmarts.keys() and bindex in atomindextoallsmarts.keys():
            asmartslist=atomindextoallsmarts[aindex]
            bsmartslist=atomindextoallsmarts[bindex]   
            combs = list(itertools.product(asmartslist, bsmartslist)) 
            for comb in combs:
                asmarts=comb[0]
                bsmarts=comb[1]
                aclass=smartstoatomclass[asmarts]
                bclass=smartstoatomclass[bsmarts]
                bondclasses=tuple([aclass,bclass])
                smartslist=[asmarts,bsmarts]
                if bondclasses not in bondclassestolistofsmartslist.keys():
                    bondclassestolistofsmartslist[bondclasses]=[]
                if tuple(bond) not in bondindicestolistofbondclasses.keys(): 
                    bondindicestolistofbondclasses[tuple(bond)]=[]
                if smartslist not in bondclassestolistofsmartslist[bondclasses]: 
                    bondclassestolistofsmartslist[bondclasses].append(smartslist)   
                 
                if bondclasses not in bondindicestolistofbondclasses[tuple(bond)]: 
                    bondindicestolistofbondclasses[tuple(bond)].append(bondclasses)    

                if bond in planarbonds or bond[::-1] in planarbonds:
                    if bondclasses not in opbendclassestolistofsmartslist.keys():
                        opbendclassestolistofsmartslist[bondclasses]=[]
                    if tuple(bond) not in opbendindicestolistofopbendclasses.keys():
                        opbendindicestolistofopbendclasses[tuple(bond)]=[]
                    if smartslist not in opbendclassestolistofsmartslist[bondclasses]: 
                        opbendclassestolistofsmartslist[bondclasses].append(smartslist)   
                    if bondclasses not in opbendindicestolistofopbendclasses[tuple(bond)]:
                        opbendindicestolistofopbendclasses[tuple(bond)].append(bondclasses)    


    for angle in listofanglesforprm:
        aindex=angle[0]
        bindex=angle[1] 
        cindex=angle[2] 
        if aindex in atomindextoallsmarts.keys() and bindex in atomindextoallsmarts.keys() and cindex in atomindextoallsmarts.keys():
            asmartslist=atomindextoallsmarts[aindex]
            bsmartslist=atomindextoallsmarts[bindex]   
            csmartslist=atomindextoallsmarts[cindex]   
            combs = list(itertools.product(asmartslist,bsmartslist,csmartslist)) 
            for comb in combs:
                asmarts=comb[0]
                bsmarts=comb[1]
                csmarts=comb[2]
                aclass=smartstoatomclass[asmarts]
                bclass=smartstoatomclass[bsmarts]
                cclass=smartstoatomclass[csmarts]
                angleclasses=tuple([aclass,bclass,cclass])
                smartslist=[asmarts,bsmarts,csmarts]
                if angleclasses not in angleclassestolistofsmartslist.keys():
                    angleclassestolistofsmartslist[angleclasses]=[]
                    strbndclassestolistofsmartslist[angleclasses]=[]
                if tuple(angle) not in angleindicestolistofangleclasses.keys():
                    angleindicestolistofangleclasses[tuple(angle)]=[]
                if tuple(angle) not in strbndindicestolistofstrbndclasses.keys():
                    strbndindicestolistofstrbndclasses[tuple(angle)]=[]
                if smartslist not in angleclassestolistofsmartslist[angleclasses]: 
                    angleclassestolistofsmartslist[angleclasses].append(smartslist)     
                if smartslist not in strbndclassestolistofsmartslist[angleclasses]:
                    strbndclassestolistofsmartslist[angleclasses].append(smartslist)     
                if angleclasses not in angleindicestolistofangleclasses[tuple(angle)]: 
                    angleindicestolistofangleclasses[tuple(angle)].append(angleclasses)     
                if angleclasses not in strbndindicestolistofstrbndclasses[tuple(angle)]: 
                    strbndindicestolistofstrbndclasses[tuple(angle)].append(angleclasses)     

    return bondclassestolistofsmartslist,angleclassestolistofsmartslist,strbndclassestolistofsmartslist,opbendclassestolistofsmartslist,bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses

def SearchForParametersViaCommentsPolarize(poltype,atomcommentstolistofsmartslist,atomindicestolistofatomcomments):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atomcommentstoparameters={}
    temp=open(poltype.latestsmallmoleculepolarizeprmlib,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if '#' not in line:
            comment=linesplit[0]
            prm=linesplit[1]
            atomcommentstoparameters[tuple([comment])]=[prm]

    atomindicestolistofatomcomments=RemoveIndicesThatDontHaveParameters(poltype,atomindicestolistofatomcomments,atomcommentstolistofsmartslist,atomcommentstoparameters,prin=True)
    return atomindicestolistofatomcomments,atomcommentstoparameters


def SearchForParametersViaCommentsNonBondedAMOEBAPlus(poltype,atomcommentstolistofsmartslist,atomindicestolistofatomcomments):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atomcommentstocpparameters={}
    atomcommentstoctparameters={}
    atomcommentstovdwparameters={}
    atomindicestolistofatomcommentscp=atomindicestolistofatomcomments.copy()
    atomindicestolistofatomcommentsct=atomindicestolistofatomcomments.copy()
    atomindicestolistofatomcommentsvdw=atomindicestolistofatomcomments.copy()

    temp=open(poltype.amoebaplusnonbondedprmlib,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        first=linesplit[0]
        if '#' not in first:
            comment=linesplit[0]
            cpprms=linesplit[1:2+1]
            linesplit[3]=str(float(linesplit[3])*3.01147) # conversion factor for ct
            ctprms=linesplit[3:4+1]
            vdwprms=linesplit[5:7+1]
            atomcommentstocpparameters[tuple([comment])]=cpprms
            atomcommentstoctparameters[tuple([comment])]=ctprms
            atomcommentstovdwparameters[tuple([comment])]=vdwprms


    atomindicestolistofatomcommentscp=RemoveIndicesThatDontHaveParameters(poltype,atomindicestolistofatomcommentscp,atomcommentstolistofsmartslist,atomcommentstocpparameters)
    atomindicestolistofatomcommentsct=RemoveIndicesThatDontHaveParameters(poltype,atomindicestolistofatomcommentsct,atomcommentstolistofsmartslist,atomcommentstoctparameters)
    atomindicestolistofatomcommentsvdw=RemoveIndicesThatDontHaveParameters(poltype,atomindicestolistofatomcommentsvdw,atomcommentstolistofsmartslist,atomcommentstovdwparameters)



    return atomindicestolistofatomcommentscp,atomcommentstocpparameters,atomindicestolistofatomcommentsct,atomcommentstoctparameters,atomindicestolistofatomcommentsvdw,atomcommentstovdwparameters


def SearchForParametersViaCommentsChargeFlux(poltype,bondcommentstolistofsmartslist,bondindicestolistofbondcomments,anglecommentstolistofsmartslist,angleindicestolistofanglecomments):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    bondcommentstocfparameters={}
    anglecommentstocfparameters={}
    temp=open(poltype.amoebapluscfprmlib,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)>0:
            first=linesplit[0]
            if '#' not in first:
                if 'cflux-b' in line:
                    comments=linesplit[2:]
                    prms=[linesplit[1]] 
                    bondcommentstocfparameters[tuple(comments)]=prms
                elif 'cflux-a' in line:
                    prms=linesplit[1:4+1]
                    comments=linesplit[5:]
                    anglecommentstocfparameters[tuple(comments)]=prms

    bondindicestolistofbondcomments=RemoveIndicesThatDontHaveParameters(poltype,bondindicestolistofbondcomments,bondcommentstolistofsmartslist,bondcommentstocfparameters)
    angleindicestolistofanglecomments=RemoveIndicesThatDontHaveParameters(poltype,angleindicestolistofanglecomments,anglecommentstolistofsmartslist,anglecommentstocfparameters)

    return bondindicestolistofbondcomments,bondcommentstocfparameters,angleindicestolistofanglecomments,anglecommentstocfparameters




def SearchForParameters(poltype,bondclassestolistofsmartslist,angleclassestolistofsmartslist,strbndclassestolistofsmartslist,opbendclassestolistofsmartslist,bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    bondclassestoparameters={}
    angleclassestoparameters={}
    strbndclassestoparameters={}
    opbendclassestoparameters={}
    temp=open(poltype.latestsmallmoleculeprmlib,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if 'bond' in line:
            bondclasses=tuple([int(linesplit[1]),int(linesplit[2])]) 
            prms=[float(linesplit[3]),float(linesplit[4])]
            bondclassestoparameters[bondclasses]=prms
        elif 'angle' in line:
            angleclasses=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]) 
            prms=[float(linesplit[4]),float(linesplit[5])]
            angleclassestoparameters[angleclasses]=prms
        elif 'strbnd' in line:
            angleclasses=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]) 
            prms=[float(linesplit[4]),float(linesplit[5])]
            strbndclassestoparameters[angleclasses]=prms
        elif 'opbend' in line:
            bondclasses=tuple([int(linesplit[1]),int(linesplit[2])]) 
            prms=[float(linesplit[5])]
            opbendclassestoparameters[bondclasses]=prms
    bondindicestolistofbondclasses=RemoveIndicesThatDontHaveParameters(poltype,bondindicestolistofbondclasses,bondclassestolistofsmartslist,bondclassestoparameters)

    angleindicestolistofangleclasses=RemoveIndicesThatDontHaveParameters(poltype,angleindicestolistofangleclasses,angleclassestolistofsmartslist,angleclassestoparameters)
    strbndindicestolistofstrbndclasses=RemoveIndicesThatDontHaveParameters(poltype,strbndindicestolistofstrbndclasses,strbndclassestolistofsmartslist,strbndclassestoparameters)
    opbendindicestolistofopbendclasses=RemoveIndicesThatDontHaveParameters(poltype,opbendindicestolistofopbendclasses,opbendclassestolistofsmartslist,opbendclassestoparameters)

    return bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses,bondclassestoparameters,angleclassestoparameters,strbndclassestoparameters,opbendclassestoparameters

def RemoveIndicesThatDontHaveParameters(poltype,indicestolistofclasses,classestolistofsmartslist,classestoparameters,prin=False):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indicestodelete=[]
    for i in range(len(indicestolistofclasses.keys())):
        indices=list(indicestolistofclasses.keys())[i]
        listofclasses=indicestolistofclasses[indices]
        classesmissing=[]
        classesexisting=[]
        
        for classes in listofclasses:
            if classes not in classestoparameters.keys() and classes[::-1] not in classestoparameters.keys():
                classesmissing.append(classes)
               
            else:
                classesexisting.append(classes)
        if len(classesmissing)==len(listofclasses):
            indicestodelete.append(indices)
        indicestolistofclasses[indices]=classesexisting
    for indices in indicestodelete:
        del indicestolistofclasses[indices]
    return indicestolistofclasses

def FindBestSMARTSMatch(poltype,indicestolistofclasses,classestolistofsmartslist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indicestoclasses={}
    indicestosmartslist={}
    for indices,listofclasses in indicestolistofclasses.items():
        fragmentsarray=[]
        classesarray=[]
        smartslistarray=[] 
        for classes in listofclasses:
            listofsmartslist=classestolistofsmartslist[classes]
            for smartslist in listofsmartslist:
                matchedindices=MatchAllSMARTS(poltype,smartslist,indices) 
                fragmentmol=CreateFragment(poltype,matchedindices)
                fragmentsarray.append(fragmentmol)
                classesarray.append(classes)
                smartslistarray.append(smartslist)
        classes,smartslist=DetermineBestSMARTSMatch(poltype,fragmentsarray,classesarray,smartslistarray)         
        indicestoclasses[indices]=classes
        indicestosmartslist[indices]=smartslist

    return indicestoclasses,indicestosmartslist

def DetermineBestSMARTSMatch(poltype,listoffragments,listofclasses,listofsmartslist):  
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tanimotoarray=[]     
    for i in range(len(listoffragments)):
        fragment=listoffragments[i]
        classes=listofclasses[i]
        smartslist=listofsmartslist[i]
        ms=[poltype.rdkitmol,fragment]
        fps = [Chem.RDKFingerprint(x) for x in ms]
        tanimoto=DataStructs.FingerprintSimilarity(fps[0],fps[1], metric=DataStructs.DiceSimilarity) 
        tanimotoarray.append(tanimoto)
    maxvalue=max(tanimotoarray)
    nums=tanimotoarray.count(maxvalue)
    if nums>1:
        smartslisttocompare=[]
        for i in range(len(listofsmartslist)):
            smartslist=listofsmartslist[i]
            tanimoto=tanimotoarray[i]
            if tanimoto==maxvalue:
                smartslisttocompare.append(smartslist)
        smartslist=ChooseMostDescriptiveSMARTSList(poltype,smartslisttocompare) 
        index=listofsmartslist.index(smartslist)
        classes=listofclasses[index]
    else:
        maxindex=tanimotoarray.index(maxvalue)
        classes=listofclasses[maxindex]
        smartslist=listofsmartslist[maxindex]
    return classes,smartslist 


def ChooseMostDescriptiveSMARTSList(poltype,smartslisttocompare):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    lengtharray=[]
    for smartslist in smartslisttocompare:
        lengths=[len(e) for e in smartslist]
        total=sum(lengths)
        lengtharray.append(total)
    maxlength=max(lengtharray)
    maxidx=lengtharray.index(maxlength)
    smartslist=smartslisttocompare[maxidx]
    return smartslist


def CreateFragment(poltype,matchedindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    m = Chem.Mol()
    em = Chem.EditableMol(m)
    oldindextonewindex={}
    for i,idx in enumerate(matchedindices):
        oldatom=poltype.rdkitmol.GetAtomWithIdx(idx)
        newindex=em.AddAtom(oldatom)
        oldindextonewindex[idx]=newindex

    for bond in poltype.rdkitmol.GetBonds():
        oendidx = bond.GetEndAtomIdx()
        obgnidx = bond.GetBeginAtomIdx()
        if oendidx not in oldindextonewindex.keys() and obgnidx not in oldindextonewindex.keys():
            continue
        elif oendidx in oldindextonewindex.keys() and obgnidx not in oldindextonewindex.keys():
            continue
        elif oendidx not in oldindextonewindex.keys() and obgnidx in oldindextonewindex.keys():
            continue

        endidx=oldindextonewindex[oendidx]
        bgnidx=oldindextonewindex[obgnidx]
        bondorder=bond.GetBondType()
        em.AddBond(bgnidx,endidx,bondorder)

    return em.GetMol()



def MatchAllSMARTS(poltype,smartslist,indices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    matchedindices=[]
    for idx in range(len(smartslist)):
        smarts=smartslist[idx]
        index=indices[idx]
        substructure = Chem.MolFromSmarts(smarts)
        diditmatch=poltype.rdkitmol.HasSubstructMatch(substructure)
        if diditmatch==True:
            matches=poltype.rdkitmol.GetSubstructMatches(substructure)
            for match in matches:
                if index in match:
                    for atomindex in match:
                        if atomindex not in matchedindices:
                            matchedindices.append(atomindex)
                    break
        else:
            raise ValueError('SMARTS match does not exist')


    return matchedindices


def GrabNewParameters(poltype,indicestoclasses,classestoparameters,keyword,indicestosmartslist,atomclasstocomment,opbendbondindicestotrigonalcenterbools=None):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indicestopoltypeclasses={}
    newpoltypeclassestocomments={}
    newpoltypeclassestosmartslist={}

    newprms=[]
    for indices,classes in indicestoclasses.items():
        smartslist=indicestosmartslist[indices]
        comments=[atomclasstocomment[k] for k in classes]
        babelindices=[k+1 for k in indices]
        if opbendbondindicestotrigonalcenterbools!=None:
            if indices in opbendbondindicestotrigonalcenterbools.keys():
                tribools=opbendbondindicestotrigonalcenterbools[indices]
            elif indices[::-1] in opbendbondindicestotrigonalcenterbools.keys():
                tribools=opbendbondindicestotrigonalcenterbools[indices[::-1]]
            if tribools[1]==False: # need second one to be trigonal center:
                babelindices=babelindices[::-1]
                indices=indices[::-1]

 
        poltypeclasses=[poltype.idxtosymclass[k] for k in babelindices]
        indicestopoltypeclasses[indices]=poltypeclasses
        newpoltypeclassestocomments[tuple([tuple(poltypeclasses)])]=comments
        newpoltypeclassestosmartslist[tuple([tuple(poltypeclasses)])]=smartslist
        poltypeclasses=[str(k) for k in poltypeclasses]
        if keyword=='opbend':
            poltypeclasses.append('0')
            poltypeclasses.append('0')
        if classes in classestoparameters.keys():
            parameters=classestoparameters[classes]
        elif classes[::-1] in classestoparameters.keys():
            parameters=classestoparameters[classes[::-1]]

        parameters=[str(k) for k in parameters]
        poltypeclassesstring=' '.join(poltypeclasses)
        parametersstring=' '.join(parameters)
        newprm=keyword+' '+poltypeclassesstring+' '+parametersstring+'\n'
        newprms.append(newprm)
        if opbendbondindicestotrigonalcenterbools!=None:
            if tribools[0]==True and tribools[1]==True:
                babelindices=babelindices[::-1]
                indices=indices[::-1]
                poltypeclasses=[poltype.idxtosymclass[k] for k in babelindices]
                indicestopoltypeclasses[indices]=poltypeclasses
                newpoltypeclassestocomments[tuple([tuple(poltypeclasses)])]=comments
                newpoltypeclassestosmartslist[tuple([tuple(poltypeclasses)])]=smartslist
                poltypeclasses=[str(k) for k in poltypeclasses]
                if keyword=='opbend':
                    poltypeclasses.append('0')
                    poltypeclasses.append('0')
                parameters=[str(k) for k in parameters]
                poltypeclassesstring=' '.join(poltypeclasses)
                parametersstring=' '.join(parameters)
                newprm=keyword+' '+poltypeclassesstring+' '+parametersstring+'\n'
                newprms.append(newprm)

        

    return indicestopoltypeclasses,newprms,newpoltypeclassestocomments,newpoltypeclassestosmartslist

def GrabNewParametersPolarize(poltype,indicestoclasses,classestoparameters,keyword,indicestosmartslist,smartstoactualcomment):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indicestopoltypeclasses={}
    newpoltypeclassestocomments={}
    newpoltypeclassestosmartslist={}

    newprms=[]
    for indices,classes in indicestoclasses.items():
        smartslist=indicestosmartslist[indices]
        comments=[smartstoactualcomment[k] for k in smartslist]
        babelindices=[k+1 for k in indices]
        poltypeclasses=[poltype.idxtosymclass[k] for k in babelindices]
        indicestopoltypeclasses[indices]=poltypeclasses
        newpoltypeclassestocomments[tuple([tuple(poltypeclasses)])]=comments
        newpoltypeclassestosmartslist[tuple([tuple(poltypeclasses)])]=smartslist
        poltypeclasses=[str(k) for k in poltypeclasses]
        if keyword=='opbend':
            poltypeclasses.append('0')
            poltypeclasses.append('0')
        if classes in classestoparameters.keys():
            parameters=classestoparameters[classes]
        elif classes[::-1] in classestoparameters.keys():
            parameters=classestoparameters[classes[::-1]]

        parameters=[str(k) for k in parameters]
        poltypeclassesstring=' '.join(poltypeclasses)
        parametersstring=' '.join(parameters)
        newprm=keyword+' '+poltypeclassesstring+' '+parametersstring+'\n'
        newprms.append(newprm)
        

    return indicestopoltypeclasses,newprms,newpoltypeclassestocomments,newpoltypeclassestosmartslist




def RemoveOldParametersKeepNewParameters(poltype,prms,newindicestopoltypeclasses,keyword,poltypeclassestoparametersmartsatomorders,poltypeclassestosmartsatomorders,poltypeclassestoelementtinkerdescrips,removefromdic=True):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newpoltypeclassestoparametersmartsatomorders=copy.deepcopy(poltypeclassestoparametersmartsatomorders)
    newpoltypeclassestosmartsatomorders=copy.deepcopy(poltypeclassestosmartsatomorders)
    newpoltypeclassestoelementtinkerdescrips=copy.deepcopy(poltypeclassestoelementtinkerdescrips) 
    newprms=[]
    dellist=[]
    appendlist=[]
    for poltypeclassesls in newpoltypeclassestoparametersmartsatomorders.keys():
        newpoltypeclasseslist=[]
        for poltypeclasses in poltypeclassesls:
            poltypeclasses=list(poltypeclasses)
            if poltypeclasses in newindicestopoltypeclasses.values() or poltypeclasses[::-1] in newindicestopoltypeclasses.values():
                pass
            else:
                newpoltypeclasseslist.append(tuple(poltypeclasses))
        newpoltypeclasseslist=tuple(newpoltypeclasseslist)
        if poltypeclassesls not in dellist:
            dellist.append(poltypeclassesls)
            if len(newpoltypeclasseslist)!=0:
                appendlist.append(newpoltypeclasseslist)


    for prm in prms:
        prmsplit=prm.split()
        if keyword=='bond' or keyword=='opbend':
            endindex=2
        elif keyword=='angle' or keyword=='strbnd':
            endindex=3
        prmsplit=prmsplit[1:endindex+1]
        classes=[int(k) for k in prmsplit]
        if classes in newindicestopoltypeclasses.values() or classes[::-1] in newindicestopoltypeclasses.values():
            pass
        else:
            newprms.append(prm)
    for i in range(len(dellist)):
        item=dellist[i]
        parametersmartatomordervalue=newpoltypeclassestoparametersmartsatomorders[item]
        smartatomordervalue=newpoltypeclassestosmartsatomorders[item]
        descrips=newpoltypeclassestoelementtinkerdescrips[item]
        if item in newpoltypeclassestoparametersmartsatomorders.keys():
            del newpoltypeclassestoparametersmartsatomorders[item]
        if item in newpoltypeclassestosmartsatomorders.keys():
            del newpoltypeclassestosmartsatomorders[item]
        if item in newpoltypeclassestoelementtinkerdescrips.keys():
            del newpoltypeclassestoelementtinkerdescrips[item]
    for newitem in appendlist:
        if len(newitem)!=0:
            newpoltypeclassestoparametersmartsatomorders[newitem]=parametersmartatomordervalue
            newpoltypeclassestosmartsatomorders[newitem]=smartatomordervalue
            newpoltypeclassestoelementtinkerdescrips[newitem]=descrips
     
    return newprms,newpoltypeclassestoparametersmartsatomorders,newpoltypeclassestosmartsatomorders,newpoltypeclassestoelementtinkerdescrips

def MapSMARTSToComments(poltype,smartstoatomclasspolar,atomclasstocommentpolar):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    smartstocomment={}
    for smarts,atomclass in smartstoatomclasspolar.items():
        comment=atomclasstocommentpolar[atomclass]
        smartstocomment[smarts]=comment

    return smartstocomment

def GrabSmartsToSoluteRadiiMap(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    smartstosoluteradiiprms={}
    temp=open(poltype.smartstosoluteradiimap,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        smarts=linesplit[0]
        prms=linesplit[1:]
        smartstosoluteradiiprms[smarts]=prms
    return smartstosoluteradiiprms


def MakeListOfListValues(poltype,keytovalue):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    newkeytovalue={}
    for key,value in keytovalue.items():
        newvalues=[]
        for item in value:
            newitem=[item]
            newvalues.append(newitem)
        newkeytovalue[key]=newvalues    

    return newkeytovalue

def GenerateSoluteParameters(poltype,atomindicestosmartslist,smartstosoluteradiiprms):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    soluteprms=[]
    for atomindices,smartslist in atomindicestosmartslist.items():
        atomindex=atomindices[0]
        smarts=smartslist[0]
        babelatomindex=atomindex+1
        poltypeclass=poltype.idxtosymclass[babelatomindex]
        prms=smartstosoluteradiiprms[smarts]
        prmsline=' '.join(prms)
        commentline='#SOLUTE-SMARTS'+' '+str(poltypeclass)+' '+smarts+'\n'
        line='SOLUTE'+' '+str(poltypeclass)+' '+prmsline+'\n'
        if line not in soluteprms:
            soluteprms.append(commentline)
            soluteprms.append(line)

    return soluteprms

def RemoveIndicesMatchedFromNewDatabase(poltype,indicestosmartsatomorders,indicestopoltypeclasses):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    removelist=[]
    for indices,smartsatomorders in indicestosmartsatomorders.items():
        if indices in indicestopoltypeclasses.keys() or indices[::-1] in indicestopoltypeclasses.keys():
            removelist.append(indices)
    for indices in removelist:
        del indicestosmartsatomorders[indices]


    return indicestosmartsatomorders


def TrimDictionary(poltype,dictotrim,dicref):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    keylist=[]
    for key,value in dictotrim.items():
        if key not in dicref.keys() and key[::-1] not in dicref.keys():
            keylist.append(key)
    for key in keylist:
        del dictotrim[key]
    return dictotrim

def AddDictionaryItems(poltype,dictoaddto,dicref):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for key,value in dicref.items():
        dictoaddto[key]=value
    return dictoaddto


def RemovePoltypeClassesFromNewMatches(poltype,missingtinkerclassestopoltypeclasses,poltypeclassestoparametersmartsatomorders):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    newmissingtinkerclassestopoltypeclasses={}
    for missingtinkerclasses,poltypeclasses in missingtinkerclassestopoltypeclasses.items():
        for polcls in poltypeclasses:
            tuppolcls=tuple([tuple(polcls)])
            if tuppolcls in poltypeclassestoparametersmartsatomorders.keys() or tuppolcls[::-1] in poltypeclassestoparametersmartsatomorders.keys():
                if missingtinkerclasses not in newmissingtinkerclassestopoltypeclasses.keys():
                    newmissingtinkerclassestopoltypeclasses[missingtinkerclasses]=[]
                newmissingtinkerclassestopoltypeclasses[missingtinkerclasses].append(polcls)
  

    return newmissingtinkerclassestopoltypeclasses

def ExtractMissingPoltypeClasses(poltype,missingtinkerclassestopoltypeclasses):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    missingpoltypeclasses=[]
    for tinkerclasses,poltypeclassesls in missingtinkerclassestopoltypeclasses.items():
        for cls in poltypeclassesls:
            missingpoltypeclasses.append(list(cls))
    return missingpoltypeclasses 


def MergeDictionaries(poltype,add,addto):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for key,value in add.items():
        addto[key]=value

    return addto


def GetPolarIndexToPolarizePrm(poltype,polarprmstotransferinfo):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    polarindextopolarizeprm={}
    for line in polarprmstotransferinfo.keys():
        linesplit=line.split()
        typenum=int(linesplit[1])
        prm=float(linesplit[2])
        for index in poltype.idxtosymclass:
            othertypenum=poltype.idxtosymclass[index]
            if typenum==othertypenum:
                polarindextopolarizeprm[index]=prm   
         

    return polarindextopolarizeprm


def TestBondAngleEquilValues(poltype):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tmpkeyfile='testbondangleequilvalues.key'
    tmpxyzfile='testbondangleequilvalues.xyz'
    alzout='testbondangleequilvaluesalz.out'
    shutil.copy(poltype.key4fname,tmpkeyfile)
    shutil.copy(poltype.xyzoutfile,tmpxyzfile)
    cmd = poltype.minimizeexe+' '+tmpxyzfile+' -k '+tmpkeyfile+' 0.1 > testbondangleequilvalues.out'
    poltype.call_subsystem([cmd], True)

    cmd = poltype.analyzeexe+' '+tmpxyzfile+'_2'+' -k '+tmpkeyfile+' '+' d > '+alzout
    poltype.call_subsystem([cmd], True)

    bondindicestonewbondequilvalues,angleindicestonewbondequilvalues=CheckBondAngleDeviationsFromQM(poltype,alzout)
   
    bondtypeindicestonewbondequilvalues=ConvertIndicesToTypeIndices(poltype,bondindicestonewbondequilvalues)
    angletypeindicestonewbondequilvalues=ConvertIndicesToTypeIndices(poltype,angleindicestonewbondequilvalues)
    ChangeBondAngleEquilValues(poltype,bondtypeindicestonewbondequilvalues,angletypeindicestonewbondequilvalues)

def ChangeBondAngleEquilValues(poltype,bondtypeindicestonewbondequilvalues,angletypeindicestonewbondequilvalues): 
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(poltype.key4fname,'r')
    results=temp.readlines()
    temp.close()
    tempname=poltype.key4fname.replace('.key','_TEMP.key')
    temp=open(tempname,'w')
    for line in results:
        if '#' not in line:
            found=False
            linesplit=line.split()
            if 'angle' in line and poltype.writeoutangle==True:
                typeindices=tuple([int(linesplit[1]),int(linesplit[2]),int(linesplit[3])])
                if typeindices in angletypeindicestonewbondequilvalues.keys():
                    equilvalue=angletypeindicestonewbondequilvalues[typeindices]
                    found=True
                elif typeindices[::-1] in angletypeindicestonewbondequilvalues.keys():
                    equilvalue=angletypeindicestonewbondequilvalues[typeindices[::-1]]
                    found=True
                if found==True:
                    linesplit=re.split(r'(\s+)', line)
                    linesplit=linesplit[:11]
                    linesplit.append('\n')
                    linesplit[10]=str(equilvalue)
                    line=''.join(linesplit)
                
            if 'bond' in line and poltype.writeoutbond==True:
                typeindices=tuple([int(linesplit[1]),int(linesplit[2])])
                if typeindices in bondtypeindicestonewbondequilvalues.keys():
                    equilvalue=bondtypeindicestonewbondequilvalues[typeindices]
                    found=True
                elif typeindices[::-1] in bondtypeindicestonewbondequilvalues.keys():
                    equilvalue=bondtypeindicestonewbondequilvalues[typeindices[::-1]]
                    found=True
                if found==True:
                    linesplit=re.split(r'(\s+)', line)
                    linesplit=linesplit[:11]
                    linesplit.append('\n')
                    linesplit[8]=str(equilvalue)
                    line=''.join(linesplit)
        temp.write(line)
    temp.close()
    os.remove(poltype.key4fname)
    os.rename(tempname,poltype.key4fname)
             


def ConvertIndicesToTypeIndices(poltype,indicestonewequilvalues):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    typeindicestonewequilvalues={}
    for indices,newequilvalues in indicestonewequilvalues.items():
        typeindices=[]
        for index in indices:
            symclass=poltype.idxtosymclass[index]
            typeindices.append(symclass)
        typeindicestonewequilvalues[tuple(typeindices)]=newequilvalues


    return typeindicestonewequilvalues


def CheckBondAngleDeviationsFromQM(poltype,alzout):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    bondindicestonewbondequilvalues={}
    angleindicestonewbondequilvalues={}
    temp=open(alzout,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if 'Additional' in line or 'Classes' in line or 'Individual' in line or 'Atom' in line or 'Stretch' in line or 'Bending' in line or 'Torsional' in line:
            continue
        linesplit=line.split()
        tol=None
        keyword=None
        if 'Angle' in line:
            indexlist=[1,2,3]
            equilindices=[4,5]
            indices,qmequil,currentequil=GrabIndicesAndEquilValues(linesplit,indexlist,equilindices)
            tol=1
            keyword='angle'
    
        if 'Bond' in line:
            indexlist=[1,2]
            equilindices=[3,4]
            indices,qmequil,currentequil=GrabIndicesAndEquilValues(linesplit,indexlist,equilindices)
            tol=2
            keyword='bond'
            
        if tol!=None:
            deviation=(np.abs(qmequil-currentequil)*100)/qmequil
            if deviation>=tol:
                shift=qmequil-(currentequil-qmequil)
                if keyword=='bond':
                    bondindicestonewbondequilvalues[tuple(indices)]=shift
                elif keyword=='angle':
                    angleindicestonewbondequilvalues[tuple(indices)]=shift



    return bondindicestonewbondequilvalues,angleindicestonewbondequilvalues

def GrabIndicesAndEquilValues(linesplit,indexlist,equilindices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    indices=[]
    for index in indexlist:
        indexele=linesplit[index] 
        indexelesplit=indexele.split('-')
        theindex=int(indexelesplit[0])
        indices.append(theindex)
    qmequil=float(linesplit[equilindices[0]])
    currentequil=float(linesplit[equilindices[1]])

    return indices,qmequil,currentequil

def RemoveIndicesMatchedFromExternalDatabase(poltype,indicestosmartsatomorder,indicestoextsmarts):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    missingindicestosmartsatomorders={}
    for indices,smartsatomorder in indicestosmartsatomorder.items():
        if indices in indicestoextsmarts.keys() or indices[::-1] in indicestoextsmarts.keys():
            pass
        else:
            missingindicestosmartsatomorders[indices]=smartsatomorder
    return missingindicestosmartsatomorders


def ExtractTransferInfo(poltype,polarprmstotransferinfo):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    polartypetotransferinfo={}
    for polarprms,transferinfo in polarprmstotransferinfo.items():
        linesplit=polarprms.split()
        polartype=linesplit[1]
        polartypetotransferinfo[polartype]=transferinfo

    return polartypetotransferinfo


def StiffenZThenBisectorAngleConstants(poltype,keyfilename):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    temp=open(keyfilename,'r')
    results=temp.readlines()
    temp.close()
    tempname=keyfilename.replace('.key','TEMP.key')
    temp=open(tempname,'w')
    angles=[]
    for line in results:
        if 'multipole ' in line and '#' not in line:
            linesplit=line.split()
            if len(linesplit)==6:
                newlinesplit=linesplit[1:-2]
                newlinesplit=[int(i) for i in newlinesplit]
                newlinesplit=[int(np.abs(i)) for i in newlinesplit]
                a=newlinesplit[0]
                b=newlinesplit[1]
                c=newlinesplit[2]
                angle=tuple([b,a,c])
                angles.append(angle)
    for line in results:
        if 'angle ' in line and '#' not in line:
            linesplit=line.split()
            indices=[linesplit[1],linesplit[2],linesplit[3]]
            indices=[int(i) for i in indices]
            indices=tuple(indices)
            for angle in angles:
                if angle==indices or angle==indices[::-1]:
                    value=float(linesplit[4])
                    value=str(2*value)
                    linesplit[4]=value
                    string='# Doubling angle force constant because z-then-bisector frame '+'\n'
                    line=string+' '.join(linesplit)
        temp.write(line)
                    
                    


    temp.close()
    os.remove(keyfilename)
    os.rename(tempname,keyfilename)


def GenerateRdkitMolObjectsParameterSMARTS(poltype,parametersmartslist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    parametersmartstordkitmol={}
    temptotalcharge=poltype.totalcharge
    for parametersmarts in parametersmartslist:
        prmmol=Chem.MolFromSmarts(parametersmarts)
        prmmol.UpdatePropertyCache()
        poltype.totalcharge=None
        prmmol,atomindextoformalcharge=poltype.CheckInputCharge(prmmol)
        Chem.SanitizeMol(prmmol)
        poltype.totalcharge=temptotalcharge
        parametersmartstordkitmol[parametersmarts]=prmmol


    return parametersmartstordkitmol


def GrabVdwParametersFromParent(poltype,oldvdwprms):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    currentdir=os.getcwd()
    parenttypestofragtypes=json.load(open("parentsymclasstofragsymclasses.txt"))
    vdwkey='parentvdw.key'
    parenttypes=list(parenttypestofragtypes.keys())
    parenttypes=[str(i) for i in parenttypes]
    fragtypes=[list(parenttypestofragtypes.values())]
    new=[]
    for fragtype in fragtypes:
        for actualtype in fragtype:
            new.append(str(actualtype[0]))
    parentvdwtransferinfo={}
    fragtypes=new[:]
    vdwprms=GrabVdwParametersFromKeyFile(poltype,vdwkey,parenttypes)
    os.chdir(currentdir)
    newvdwprms=[]
    for oldvdwprmline in oldvdwprms:
        linesplit=oldvdwprmline.split()
        typenum=linesplit[1]
        if typenum not in fragtypes:
            newvdwprms.append(oldvdwprmline)
    for vdwprmline in vdwprms:
        linesplit=vdwprmline.split()
        typenum=linesplit[1]
        fragtype=str(parenttypestofragtypes[typenum][0])
        linesplit[1]=fragtype
        otherline='# Transferring from parent vdw type '+typenum+'\n'
        parentvdwtransferinfo[fragtype]=otherline
        newline=' '.join(linesplit)+'\n'
        newvdwprms.append(newline)

    return newvdwprms,parentvdwtransferinfo


def GrabVdwParametersFromKeyFile(poltype,vdwkey,parenttypes):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    vdwprms=[]
    temp=open(vdwkey,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if 'vdw ' in line:
            typenum=linesplit[1]
            for parenttype in parenttypes:
                if str(parenttype)==typenum:
                    vdwprms.append(line)

    return vdwprms



def AddParentVdwTransferInfo(poltype,parentvdwtransferinfo,vdwprmstotransferinfo):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    for line in vdwprmstotransferinfo.keys():
        linesplit=line.split()
        typenum=linesplit[1]
        if str(typenum) in parentvdwtransferinfo.keys():
            otherline=parentvdwtransferinfo[typenum]
            vdwprmstotransferinfo[line]+=otherline



    return vdwprmstotransferinfo

def AddTrivialOPBendForAmine(poltype,opbendprms,mol):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    atomiter=openbabel.OBMolAtomIter(mol)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        atomicnum=atom.GetAtomicNum()
        if atomicnum==7:
            hcount=0
            neighbs=[natom for natom in openbabel.OBAtomAtomIter(atom)]
            atomidxtoatomicnum={}
            for neighb in neighbs:
                natomidx=neighb.GetIdx()
                natomicnum=neighb.GetAtomicNum()
                atomidxtoatomicnum[natomidx]=natomicnum
                if natomicnum==1:
                    hcount+=1
            if hcount==2:
                for natomidx,atomicnum in atomidxtoatomicnum.items():
                    if atomicnum!=1:
                        theidx=natomidx
                    elif atomicnum==1:
                        oidx=natomidx
                indices=[theidx,atomidx]
                symclasses=[poltype.idxtosymclass[i] for i in indices]
                symclasses=[str(i) for i in symclasses]
                string=' '.join(symclasses)
                string='opbend '+string+' 0 0 0'+'\n'
                opbendprms.append(string)
                indices=[oidx,atomidx]
                symclasses=[poltype.idxtosymclass[i] for i in indices]
                symclasses=[str(i) for i in symclasses]
                string=' '.join(symclasses)
                string='opbend '+string+' 0 0 0'+'\n'
                opbendprms.append(string)



    return opbendprms

def GenerateTinkerClassesToPoltypeClasses(poltype,indices):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """
    tinkerclassestopoltypeclasses={}
    for indexls in indices:
        indexls=[i+1 for i in indexls]
        poltypeclasses=[poltype.idxtosymclass[i] for i in indexls]
        tup=tuple(poltypeclasses)
        ls=[tup]
        tinkerclassestopoltypeclasses[tup]=ls

    return tinkerclassestopoltypeclasses



def ForceSameResonanceTypesSameMatches(poltype,atomindicestocomments,atomindicestosmartslist):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    for atomindices,comments in atomindicestocomments.items():
        smartslist=atomindicestosmartslist[atomindices]
        atomindex=atomindices[0]
        atomtype=poltype.idxtosymclass[atomindex+1]
        equivalentindices=[]
        for idx,typenum in poltype.idxtosymclass.items():
            if typenum==atomtype:
                equivalentindices.append(idx-1)
        equivalentindextosmartslist={}
        for idx in equivalentindices:
            oatomindices=tuple([idx])
            osmartslist=atomindicestosmartslist[oatomindices]
            equivalentindextosmartslist[oatomindices]=osmartslist
        lentoidx={}
        for eidx,smartslist in equivalentindextosmartslist.items():
            length=len(smartslist[0])
            lentoidx[length]=eidx
        maxlen=max(lentoidx.keys())
        maxidx=lentoidx[maxlen]
        maxsmartslist=atomindicestosmartslist[maxidx]
        maxcomments=atomindicestocomments[maxidx]
        atomindicestocomments[atomindices]=maxcomments
        atomindicestosmartslist[atomindices]=maxsmartslist


    return atomindicestocomments,atomindicestosmartslist


def GrabSmallMoleculeAMOEBAParameters(poltype,optmol,mol,rdkitmol,polarize=False):
    """
    Intent:
    Input:
    Output:
    Referenced By: 
    Description: 
    """

    if polarize==True:
        listofatomsforprm,listofbondsforprm,listofanglesforprm,listoftorsionsforprm=GrabAtomsForParameters(poltype,mol)

        smartstoatomclasspolar, atomclasstocommentpolar,smartstoactualcomment=ReadDatabaseSmartsMapPolarize(poltype,poltype.latestsmallmoleculesmartstotypespolarize) 
        atomclasstoatomclass=dict(zip(atomclasstocommentpolar.keys(),atomclasstocommentpolar.keys()))
       
        atomindextoallsmartspolar,atomindextoallsmartsmatchespolar=MatchAllSmartsToAtomIndices(poltype,smartstoatomclasspolar)
        smartstocomment=MapSMARTSToComments(poltype,smartstoatomclasspolar,atomclasstoatomclass)
        
        atomcommentstolistofsmartslist,atomindicestolistofatomcomments=MapIndicesToCommentsAtom(poltype,atomindextoallsmartspolar,smartstocomment,listofatomsforprm)
        atomindicestolistofatomcomments,atomcommentstoparameters=SearchForParametersViaCommentsPolarize(poltype,atomcommentstolistofsmartslist,atomindicestolistofatomcomments)
        atomindicestocomments,atomindicestosmartslist=FindBestSMARTSMatch(poltype,atomindicestolistofatomcomments,atomcommentstolistofsmartslist)

        atomindicestocomments,atomindicestosmartslist=ForceSameResonanceTypesSameMatches(poltype,atomindicestocomments,atomindicestosmartslist)




        newpolarindicestopoltypeclasses,newpolarprms,newpolarpoltypecommentstocomments,newpolarpoltypecommentstosmartslist=GrabNewParametersPolarize(poltype,atomindicestocomments,atomcommentstoparameters,'polarize',atomindicestosmartslist,smartstoactualcomment) 
        
        polarprmstotransferinfo=MapParameterLineToTransferInfo(poltype,newpolarprms,{},{},{},{},newpolarpoltypecommentstocomments,newpolarpoltypecommentstosmartslist,{},[],[],[],[])
        polarindextopolarizeprm=GetPolarIndexToPolarizePrm(poltype,polarprmstotransferinfo)
        polartypetotransferinfo=ExtractTransferInfo(poltype,polarprmstotransferinfo)
        return polarindextopolarizeprm,polartypetotransferinfo

    else:

    
        listofatomsforprm,listofbondsforprm,listofanglesforprm,listoftorsionsforprm=GrabAtomsForParameters(poltype,mol)
        if poltype.forcefield=='AMOEBA+':
            smartstoatomclasscf, atomclasstoclassnamecf, atomclasstocommentcf=ReadDatabaseSmartsMap(poltype,poltype.amoebapluscfsmartstocommentmap) 
            atomindextoallsmartscf,atomindextoallsmartsmatchescf=MatchAllSmartsToAtomIndices(poltype,smartstoatomclasscf)
            smartstocommentcf=MapSMARTSToComments(poltype,smartstoatomclasscf,atomclasstoclassnamecf)
            bondcommentstolistofsmartslist,bondindicestolistofbondcomments,anglecommentstolistofsmartslist,angleindicestolistofanglecomments=MapIndicesToCommentsBondAngle(poltype,atomindextoallsmartscf,smartstocommentcf,listofbondsforprm,listofanglesforprm)
            bondindicestolistofbondcomments,bondcommentstocfparameters,angleindicestolistofanglecomments,anglecommentstocfparameters=SearchForParametersViaCommentsChargeFlux(poltype,bondcommentstolistofsmartslist,bondindicestolistofbondcomments,anglecommentstolistofsmartslist,angleindicestolistofanglecomments)

            bondindicestocommentscf,bondindicestosmartslistcf=FindBestSMARTSMatch(poltype,bondindicestolistofbondcomments,bondcommentstolistofsmartslist)
            angleindicestocommentscf,angleindicestosmartslistcf=FindBestSMARTSMatch(poltype,angleindicestolistofanglecomments,anglecommentstolistofsmartslist)
            commentlistcf=list(atomclasstoclassnamecf.values())
            atomcommentstocommentcf=dict(zip(commentlistcf,commentlistcf))
            bondcfindicestopoltypeclasses,bondcfprms,bondcfpoltypecommentstocomments,bondcfpoltypecommentstosmartslist=GrabNewParameters(poltype,bondindicestocommentscf,bondcommentstocfparameters,'cflux-b',bondindicestosmartslistcf,atomcommentstocommentcf) 
            anglecfindicestopoltypeclasses,anglecfprms,anglecfpoltypecommentstocomments,anglecfpoltypecommentstosmartslist=GrabNewParameters(poltype,angleindicestocommentscf,anglecommentstocfparameters,'cflux-a',angleindicestosmartslistcf,atomcommentstocommentcf) 
            smartstoatomclassnonbondedplus, atomclasstoclassnamenonbondedplus, atomclasstocommentnonbondedplus=ReadDatabaseSmartsMap(poltype,poltype.amoebaplusnonbondedsmartstocommentmap) 
            atomindextoallsmartsnonbondedplus,atomindextoallsmartsmatchesnonbondedplus=MatchAllSmartsToAtomIndices(poltype,smartstoatomclassnonbondedplus)
            smartstocommentnonbondedplus=MapSMARTSToComments(poltype,smartstoatomclassnonbondedplus,atomclasstoclassnamenonbondedplus)
            atomcommentstolistofsmartslistnonbondedplus,atomindicestolistofatomcommentsnonbondedplus=MapIndicesToCommentsAtom(poltype,atomindextoallsmartsnonbondedplus,smartstocommentnonbondedplus,listofatomsforprm)

            atomindicestolistofatomcommentscp,atomcommentstocpparameters,atomindicestolistofatomcommentsct,atomcommentstoctparameters,atomindicestolistofatomcommentsvdw,atomcommentstovdwparameters=SearchForParametersViaCommentsNonBondedAMOEBAPlus(poltype,atomcommentstolistofsmartslistnonbondedplus,atomindicestolistofatomcommentsnonbondedplus)
            atomindicestocommentscp,atomindicestosmartslistcp=FindBestSMARTSMatch(poltype,atomindicestolistofatomcommentscp,atomcommentstolistofsmartslistnonbondedplus)
            atomindicestocommentsct,atomindicestosmartslistct=FindBestSMARTSMatch(poltype,atomindicestolistofatomcommentsct,atomcommentstolistofsmartslistnonbondedplus)
            atomindicestocommentsvdw,atomindicestosmartslistvdw=FindBestSMARTSMatch(poltype,atomindicestolistofatomcommentsvdw,atomcommentstolistofsmartslistnonbondedplus)

            commentlistnonbondedplus=list(atomclasstoclassnamenonbondedplus.values())
            atomcommentstocommentnonbondedplus=dict(zip(commentlistnonbondedplus,commentlistnonbondedplus))

            cpindicestopoltypeclasses,cpprms,cppoltypecommentstocomments,cppoltypecommentstosmartslist=GrabNewParameters(poltype,atomindicestocommentscp,atomcommentstocpparameters,'chgpen',atomindicestosmartslistcp,atomcommentstocommentnonbondedplus) 
            ctindicestopoltypeclasses,ctprms,ctpoltypecommentstocomments,ctpoltypecommentstosmartslist=GrabNewParameters(poltype,atomindicestocommentsct,atomcommentstoctparameters,'chgtrn',atomindicestosmartslistct,atomcommentstocommentnonbondedplus) 
            newvdwindicestopoltypeclasses,newvdwprms,newvdwpoltypecommentstocomments,newvdwpoltypecommentstosmartslist=GrabNewParameters(poltype,atomindicestocommentsvdw,atomcommentstovdwparameters,'vdw',atomindicestosmartslistvdw,atomcommentstocommentnonbondedplus) 
        smartstosoluteradiiprms=GrabSmartsToSoluteRadiiMap(poltype)   
        atomindextoallsmartssolute,atomindextoallsmartsmatchessolute=MatchAllSmartsToAtomIndices(poltype,smartstosoluteradiiprms)
        atomindices=list(atomindextoallsmartssolute.keys())
        atomindices=[tuple([k]) for k in atomindices]
        atomindicestolistofatomindices=dict(zip(atomindices,atomindices))
        atomindextoallsmartssolute=MakeListOfListValues(poltype,atomindextoallsmartssolute)
        atomindicestoindices,atomindicestosmartslist=FindBestSMARTSMatch(poltype,atomindicestolistofatomindices,atomindextoallsmartssolute)
        soluteprms=GenerateSoluteParameters(poltype,atomindicestosmartslist,smartstosoluteradiiprms)
        smartstoatomclasspolar, atomclasstocommentpolar,smartstoactualcomment=ReadDatabaseSmartsMapPolarize(poltype,poltype.latestsmallmoleculesmartstotypespolarize) 
        atomclasstoatomclass=dict(zip(atomclasstocommentpolar.keys(),atomclasstocommentpolar.keys()))
        atomindextoallsmartspolar,atomindextoallsmartsmatchespolar=MatchAllSmartsToAtomIndices(poltype,smartstoatomclasspolar)
        smartstocomment=MapSMARTSToComments(poltype,smartstoatomclasspolar,atomclasstoatomclass)
        atomcommentstolistofsmartslist,atomindicestolistofatomcomments=MapIndicesToCommentsAtom(poltype,atomindextoallsmartspolar,smartstocomment,listofatomsforprm)
        atomindicestolistofatomcomments,atomcommentstoparameters=SearchForParametersViaCommentsPolarize(poltype,atomcommentstolistofsmartslist,atomindicestolistofatomcomments)
        atomindicestocomments,atomindicestosmartslist=FindBestSMARTSMatch(poltype,atomindicestolistofatomcomments,atomcommentstolistofsmartslist)
        newpolarindicestopoltypeclasses,newpolarprms,newpolarpoltypecommentstocomments,newpolarpoltypecommentstosmartslist=GrabNewParametersPolarize(poltype,atomindicestocomments,atomcommentstoparameters,'polarize',atomindicestosmartslist,smartstoactualcomment) 
        smartstoatomclass, atomclasstoclassname, atomclasstocomment=ReadDatabaseSmartsMap(poltype,poltype.latestsmallmoleculesmartstotinkerclass) 
        atomindextoallsmarts,atomindextoallsmartsmatches=MatchAllSmartsToAtomIndices(poltype,smartstoatomclass)
        bondsmartsatomordertoparameters,anglesmartsatomordertoparameters,strbndsmartsatomordertoparameters,torsionsmartsatomordertoparameters,opbendsmartsatomordertoparameters,vdwsmartsatomordertoparameters,tortorsmartsatomordertoparameters,tortorsmartsatomordertogrid,smartsatomordertotorvdwdb=ReadExternalDatabase(poltype)
        smartsatomordertoelementtinkerdescrip=ReadSmallMoleculeLib(poltype,poltype.smallmoleculesmartstotinkerdescrip)
        elementtinkerdescriptotinkertype,tinkertypetoclass=GrabTypeAndClassNumbers(poltype,poltype.smallmoleculeprmlib)
        planarbonds=GrabPlanarBonds(poltype,listofbondsforprm,mol)

        bondclassestolistofsmartslist,angleclassestolistofsmartslist,strbndclassestolistofsmartslist,opbendclassestolistofsmartslist,bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses=MapIndicesToClasses(poltype,atomindextoallsmarts,smartstoatomclass,listofbondsforprm,listofanglesforprm,planarbonds)
        bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses,bondclassestoparameters,angleclassestoparameters,strbndclassestoparameters,opbendclassestoparameters=SearchForParameters(poltype,bondclassestolistofsmartslist,angleclassestolistofsmartslist,strbndclassestolistofsmartslist,opbendclassestolistofsmartslist,bondindicestolistofbondclasses,angleindicestolistofangleclasses,strbndindicestolistofstrbndclasses,opbendindicestolistofopbendclasses)
        bondindicestoclasses,bondindicestosmartslist=FindBestSMARTSMatch(poltype,bondindicestolistofbondclasses,bondclassestolistofsmartslist)

 
        angleindicestoclasses,angleindicestosmartslist=FindBestSMARTSMatch(poltype,angleindicestolistofangleclasses,angleclassestolistofsmartslist)
 
        strbndindicestoclasses,strbndindicestosmartslist=FindBestSMARTSMatch(poltype,strbndindicestolistofstrbndclasses,strbndclassestolistofsmartslist)
        opbendindicestoclasses,opbendindicestosmartslist=FindBestSMARTSMatch(poltype,opbendindicestolistofopbendclasses,opbendclassestolistofsmartslist)
        opbendbondindicestotrigonalcenterbools=CheckTrigonalCenters(poltype,listofbondsforprm,mol)
        newangleindicestopoltypeclasses,newangleprms,newanglepoltypeclassestocomments,newanglepoltypeclassestosmartslist=GrabNewParameters(poltype,angleindicestoclasses,angleclassestoparameters,'angle',angleindicestosmartslist,atomclasstocomment) 
        newbondindicestopoltypeclasses,newbondprms,newbondpoltypeclassestocomments,newbondpoltypeclassestosmartslist=GrabNewParameters(poltype,bondindicestoclasses,bondclassestoparameters,'bond',bondindicestosmartslist,atomclasstocomment) 
        newstrbndindicestopoltypeclasses,newstrbndprms,newstrbndpoltypeclassestocomments,newstrbndpoltypeclassestosmartslist=GrabNewParameters(poltype,strbndindicestoclasses,strbndclassestoparameters,'strbnd',strbndindicestosmartslist,atomclasstocomment) 
        newopbendindicestopoltypeclasses,newopbendprms,newopbendpoltypeclassestocomments,newopbendpoltypeclassestosmartslist=GrabNewParameters(poltype,opbendindicestoclasses,opbendclassestoparameters,'opbend',opbendindicestosmartslist,atomclasstocomment,opbendbondindicestotrigonalcenterbools) 
 
        
        listofanglesthatneedplanarkeyword=CheckForPlanerAngles(poltype,listofanglesforprm,mol)
        parametersmartslist=GrabSMARTSList(poltype,smartsatomordertoelementtinkerdescrip)
        parametersmartstomaxcommonsubstructuresmarts,maxatomsize=FindMaximumCommonSubstructures(poltype,parametersmartslist,rdkitmol)
        if len(list(parametersmartstomaxcommonsubstructuresmarts.keys()))!=0:
             parametersmartslist=list(parametersmartstomaxcommonsubstructuresmarts.keys())
        indextoneighbidxs=FindAllNeighborIndexes(poltype,rdkitmol)
        bondindicestoextsmartsmatchlength,bondindicestoextsmarts,bondindicestoextsmartsatomorder,bondsmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,bondsmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb)
        angleindicestoextsmartsmatchlength,angleindicestoextsmarts,angleindicestoextsmartsatomorder,anglesmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,anglesmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb)
        strbndindicestoextsmartsmatchlength,strbndindicestoextsmarts,strbndindicestoextsmartsatomorder,strbndsmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,strbndsmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb)
        torsionindicestoextsmartsmatchlength,torsionindicestoextsmarts,torsionindicestoextsmartsatomorders,torsionsmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,torsionsmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb)
        opbendindicestoextsmartsmatchlength,opbendindicestoextsmarts,opbendindicestoextsmartsatomorder,opbendsmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,opbendsmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb)
        vdwindicestoextsmartsmatchlength,vdwindicestoextsmarts,vdwindicestoextsmartsatomorder,vdwsmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,vdwsmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb,torsionindicestoextsmarts)
        tortorindicestoextsmartsmatchlength,tortorindicestoextsmarts,tortorindicestoextsmartsatomorders,tortorsmartsatomordertoparameters=MatchExternalSMARTSToMolecule(poltype,rdkitmol,tortorsmartsatomordertoparameters,indextoneighbidxs,smartsatomordertotorvdwdb)
        parametersmartstordkitmol=GenerateRdkitMolObjectsParameterSMARTS(poltype,parametersmartslist)
        atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,atomindicesforprmtomatchallneighbs=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,listofatomsforprm,parametersmartslist,mol,parametersmartstordkitmol)
        bondsforprmtoparametersmarts,bondsforprmtosmarts,bondsforprmtomatchallneighbs=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,listofbondsforprm,parametersmartslist,mol,parametersmartstordkitmol)
        planarbondsforprmtoparametersmarts,planarbondsforprmtosmarts,planarbondsforprmtomatchallneighbs=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,planarbonds,parametersmartslist,mol,parametersmartstordkitmol)
        anglesforprmtoparametersmarts,anglesforprmtosmarts,anglesforprmtomatchallneighbs=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,listofanglesforprm,parametersmartslist,mol,parametersmartstordkitmol)
        planaranglesforprmtoparametersmarts,planaranglesforprmtosmarts,planaranglesforprmtomatchallneighbs=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,listofanglesforprm,parametersmartslist,mol,parametersmartstordkitmol)
        torsionsforprmtoparametersmarts,torsionsforprmtosmarts,torsionsforprmtomatchallneighbs=MatchAtomIndicesSMARTSToParameterSMARTS(poltype,listoftorsionsforprm,parametersmartslist,mol,parametersmartstordkitmol)
        atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,vdwindicestoextsmarts,vdwindicestoextsmartsatomorder=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,vdwindicestoextsmartsmatchlength,atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,vdwindicestoextsmarts,atomindicesforprmtomatchallneighbs,vdwindicestoextsmartsatomorder)
        bondsforprmtoparametersmarts,bondsforprmtosmarts,bondindicestoextsmarts,bondindicestoextsmartsatomorder=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,bondindicestoextsmartsmatchlength,bondsforprmtoparametersmarts,bondsforprmtosmarts,bondindicestoextsmarts,bondsforprmtomatchallneighbs,bondindicestoextsmartsatomorder)
        planarbondsforprmtoparametersmarts,planarbondsforprmtosmarts,planarbondindicestoextsmarts,planarbondindicestoextsmartsatomorder=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,bondindicestoextsmartsmatchlength,planarbondsforprmtoparametersmarts,planarbondsforprmtosmarts,bondindicestoextsmarts,planarbondsforprmtomatchallneighbs,bondindicestoextsmartsatomorder)
        anglesforprmtoparametersmarts,anglesforprmtosmarts,angleindicestoextsmarts,angleindicestoextsmartsatomorder=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,angleindicestoextsmartsmatchlength,anglesforprmtoparametersmarts,anglesforprmtosmarts,angleindicestoextsmarts,anglesforprmtomatchallneighbs,angleindicestoextsmartsatomorder)
        planaranglesforprmtoparametersmarts,planaranglesforprmtosmarts,planarangleindicestoextsmarts,planarangleindicestoextsmartsatomorder=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,angleindicestoextsmartsmatchlength,planaranglesforprmtoparametersmarts,planaranglesforprmtosmarts,angleindicestoextsmarts,planaranglesforprmtomatchallneighbs,angleindicestoextsmartsatomorder)
        torsionsforprmtoparametersmarts,torsionsforprmtosmarts,torsionindicestoextsmarts,torsionindicestoextsmartsatomorder=CompareParameterSMARTSMatchToExternalSMARTSMatch(poltype,torsionindicestoextsmartsmatchlength,torsionsforprmtoparametersmarts,torsionsforprmtosmarts,torsionindicestoextsmarts,torsionsforprmtomatchallneighbs,torsionindicestoextsmartsatomorders)
        torsionindicestoextsmartsatomorders=TrimDictionary(poltype,torsionindicestoextsmartsatomorders,torsionindicestoextsmarts)
        atomindextotinkertype,atomindextotinkerclass,atomindextoparametersmartsatomorder,atomindextoelementtinkerdescrip,atomindextosmartsatomorder=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,atomindicesforprmtoparametersmarts,atomindicesforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)
 
        bondindicestotinkertypes,bondindicestotinkerclasses,bondindicestoparametersmartsatomorders,bondindicestoelementtinkerdescrips,bondindicestosmartsatomorders=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,bondsforprmtoparametersmarts,bondsforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)

 
        opbendbondindicestotinkerclasses,opbendbondindicestosmartsatomorders=FilterBondSMARTSEnviorment(poltype,bondindicestosmartsatomorders,bondindicestotinkerclasses)

        newplanarbonds=ConvertListOfListToListOfTuples(poltype,planarbonds)
        newdics=FilterDictionaries(poltype,[opbendbondindicestotinkerclasses,opbendbondindicestosmartsatomorders],newplanarbonds)
        opbendbondindicestotinkerclasses,opbendbondindicestosmartsatomorders=newdics[:]
        planarbondindicestotinkertypes,planarbondindicestotinkerclasses,planarbondindicestoparametersmartsatomorders,planarbondindicestoelementtinkerdescrips,planarbondindicestosmartsatomorders=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,planarbondsforprmtoparametersmarts,planarbondsforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)

        angleindicestotinkertypes,angleindicestotinkerclasses,angleindicestoparametersmartsatomorders,angleindicestoelementtinkerdescrips,angleindicestosmartsatomorders=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,anglesforprmtoparametersmarts,anglesforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)
        planarangleindicestotinkertypes,planarangleindicestotinkerclasses,planarangleindicestoparametersmartsatomorders,planarangleindicestoelementtinkerdescrips,planarangleindicestosmartsatomorders=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,planaranglesforprmtoparametersmarts,planaranglesforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)
        torsionindicestotinkertypes,torsionindicestotinkerclasses,torsionindicestoparametersmartsatomorders,torsionindicestoelementtinkerdescrips,torsionindicestosmartsatomorders=GenerateAtomIndexToAtomTypeAndClassForAtomList(poltype,torsionsforprmtoparametersmarts,torsionsforprmtosmarts,smartsatomordertoelementtinkerdescrip,elementtinkerdescriptotinkertype,tinkertypetoclass,rdkitmol)
        originalbondindicestosmartsatomorders=bondindicestosmartsatomorders.copy()
        originalangleindicestosmartsatomorders=angleindicestosmartsatomorders.copy()
        formissingangleindicestosmartsatomorders=RemoveIndicesMatchedFromNewDatabase(poltype,angleindicestosmartsatomorders,newangleindicestopoltypeclasses) 
        
        formissingbondindicestosmartsatomorders=RemoveIndicesMatchedFromNewDatabase(poltype,bondindicestosmartsatomorders,newbondindicestopoltypeclasses)
        formissingvdwindicestosmartsatomorders=RemoveIndicesMatchedFromExternalDatabase(poltype,atomindextosmartsatomorder,vdwindicestoextsmarts) 
        bondindicestosmartsatomorders=originalbondindicestosmartsatomorders.copy()
        angleindicestosmartsatomorders=originalangleindicestosmartsatomorders.copy()

        totalbondscollector=FindAllConsecutiveRotatableBonds(poltype,mol,listofbondsforprm)
        tortorsmissing=FindMissingTorTors(poltype,tortorindicestoextsmarts,tortorsmartsatomordertoparameters,rdkitmol,mol,indextoneighbidxs,totalbondscollector)
        torsionindicestosmartsatomorders=AddDictionaryItems(poltype,torsionindicestosmartsatomorders,torsionindicestoextsmartsatomorders)
        torsionsmissing,poormatchingaromatictorsions,poormatchingpartialaromatictorsions,torsionstozerooutduetocolinear=FindMissingTorsions(poltype,torsionindicestosmartsatomorders,rdkitmol,optmol,indextoneighbidxs)
        torsionsmissing=FindAdjacentMissingTorsionsForTorTor(poltype,torsionsmissing,totalbondscollector,tortorsmissing)
        atomindextosmartsatomorder=AddExternalDatabaseMatches(poltype, atomindextosmartsatomorder,vdwindicestoextsmarts,vdwsmartsatomordertoparameters)
        vdwmissing=FindMissingParameters(poltype,formissingvdwindicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs)
        missingvdwatomindices=ReduceMissingVdwByTypes(poltype,vdwmissing)
        bondmissing=FindMissingParameters(poltype,formissingbondindicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs)
        anglemissing=FindMissingParameters(poltype,formissingangleindicestosmartsatomorders,rdkitmol,mol,indextoneighbidxs)
        anglemissingindicestotinkerclasses,removedangleindices=PruneDictionary(poltype,anglemissing,angleindicestotinkerclasses)
        bondmissingindicestotinkerclasses,removedbondindices=PruneDictionary(poltype,bondmissing,bondindicestotinkerclasses)
        torsionsmissingindicestotinkerclasses,removedtorsionindices=PruneDictionary(poltype,torsionsmissing,torsionindicestotinkerclasses)
        torsionszerooutindicestotinkerclasses,removedzerotorsionindices=PruneDictionary(poltype,torsionstozerooutduetocolinear,torsionindicestotinkerclasses)
        atomtinkerclasstopoltypeclass=TinkerClassesToPoltypeClasses(poltype,atomindextotinkerclass)
        bondtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,bondindicestotinkerclasses)
        arotorsionsmissingindicestotinkerclasses,removedaromissingtorsionindices=PruneDictionary(poltype,poormatchingaromatictorsions,torsionindicestotinkerclasses)
        partialarotorsionsmissingindicestotinkerclasses,removedpartialarotorsionmpartialarotorsionmissinggindices=PruneDictionary(poltype,poormatchingpartialaromatictorsions,torsionindicestotinkerclasses)
        planarbondtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,planarbondindicestotinkerclasses,False)
        opbendtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,opbendbondindicestotinkerclasses)
        opbendtinkerclassestopoltypeclasses=AddReverseKeys(poltype,opbendtinkerclassestopoltypeclasses)
        opbendtinkerclassestotrigonalcenterbools=TinkerClassesToTrigonalCenter(poltype,opbendbondindicestotinkerclasses,opbendbondindicestotrigonalcenterbools)

        angletinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,angleindicestotinkerclasses)
        anglepoltypeclassestotinkerclasses=ReverseDictionaryValueList(poltype,angletinkerclassestopoltypeclasses)

        planarangletinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,planarangleindicestotinkerclasses)
        torsiontinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,torsionindicestotinkerclasses)
        torsionsmissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,torsionsmissingindicestotinkerclasses)
        torsionszeroouttinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,torsionszerooutindicestotinkerclasses)

        arotorsionsmissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,arotorsionsmissingindicestotinkerclasses)
        partialarotorsionsmissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,partialarotorsionsmissingindicestotinkerclasses)
        anglemissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,anglemissingindicestotinkerclasses)
        bondmissingtinkerclassestopoltypeclasses=TinkerClassesToPoltypeClasses(poltype,bondmissingindicestotinkerclasses)
        bondpoltypeclassestotinkerclasses=ReverseDictionaryValueList(poltype,bondtinkerclassestopoltypeclasses)
        torsionpoltypeclassestoparametersmartsatomorders=ConvertIndicesDictionaryToPoltypeClasses(poltype,torsionindicestoparametersmartsatomorders,torsionindicestotinkerclasses,torsiontinkerclassestopoltypeclasses)
        torsionpoltypeclassestosmartsatomorders=ConvertIndicesDictionaryToPoltypeClasses(poltype,torsionindicestosmartsatomorders,torsionindicestotinkerclasses,torsiontinkerclassestopoltypeclasses)
        torsionpoltypeclassestoelementtinkerdescrips=ConvertIndicesDictionaryToPoltypeClasses(poltype,torsionindicestoelementtinkerdescrips,torsionindicestotinkerclasses,torsiontinkerclassestopoltypeclasses)
        anglepoltypeclassestoparametersmartsatomorders=ConvertIndicesDictionaryToPoltypeClasses(poltype,angleindicestoparametersmartsatomorders,angleindicestotinkerclasses,angletinkerclassestopoltypeclasses)
        anglepoltypeclassestosmartsatomorders=ConvertIndicesDictionaryToPoltypeClasses(poltype,angleindicestosmartsatomorders,angleindicestotinkerclasses,angletinkerclassestopoltypeclasses)
        anglepoltypeclassestoelementtinkerdescrips=ConvertIndicesDictionaryToPoltypeClasses(poltype,angleindicestoelementtinkerdescrips,angleindicestotinkerclasses,angletinkerclassestopoltypeclasses)

        bondpoltypeclassestoparametersmartsatomorders=ConvertIndicesDictionaryToPoltypeClasses(poltype,bondindicestoparametersmartsatomorders,bondindicestotinkerclasses,bondtinkerclassestopoltypeclasses)
        bondpoltypeclassestosmartsatomorders=ConvertIndicesDictionaryToPoltypeClasses(poltype,bondindicestosmartsatomorders,bondindicestotinkerclasses,bondtinkerclassestopoltypeclasses)

        bondpoltypeclassestoelementtinkerdescrips=ConvertIndicesDictionaryToPoltypeClasses(poltype,bondindicestoelementtinkerdescrips,bondindicestotinkerclasses,bondtinkerclassestopoltypeclasses)
        atompoltypeclasstoparametersmartsatomorder=ConvertIndicesDictionaryToPoltypeClasses(poltype,atomindextoparametersmartsatomorder,atomindextotinkerclass,atomtinkerclasstopoltypeclass)

        atompoltypeclasstosmartsatomorder=ConvertIndicesDictionaryToPoltypeClasses(poltype,atomindextosmartsatomorder,atomindextotinkerclass,atomtinkerclasstopoltypeclass)
        atompoltypeclassestoelementtinkerdescrip=ConvertIndicesDictionaryToPoltypeClasses(poltype,atomindextoelementtinkerdescrip,atomindextotinkerclass,atomtinkerclasstopoltypeclass)
        poltypetoprmtype={} # dont need polarize parameters 
        typestoframedefforprmfile={} # dont need multipole parameters
        fname=poltype.smallmoleculeprmlib
        bondprms,angleprms,torsionprms,strbndprms,mpoleprms,opbendprms,polarizeprms,vdwprms,torsiontopitor=GrabParametersFromPrmFile(poltype,bondtinkerclassestopoltypeclasses,opbendtinkerclassestopoltypeclasses,opbendtinkerclassestotrigonalcenterbools,angletinkerclassestopoltypeclasses,torsiontinkerclassestopoltypeclasses,poltypetoprmtype,atomtinkerclasstopoltypeclass,typestoframedefforprmfile,fname,True)
        torsionprms=CorrectPitorEnergy(poltype,torsionprms,torsiontopitor)
        angleprms,anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips=RemoveOldParametersKeepNewParameters(poltype,angleprms,newangleindicestopoltypeclasses,'angle',anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips)
        opbendprms,blankbondpoltypeclassestoparametersmartsatomorders,blankbondpoltypeclassestosmartsatomorders,blankbondpoltypeclassestoelementtinkerdescrips=RemoveOldParametersKeepNewParameters(poltype,opbendprms,newopbendindicestopoltypeclasses,'opbend',bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips)
        opbendpoltypeclassestosmartsatomorders=bondpoltypeclassestosmartsatomorders.copy()
        opbendpoltypeclassestoelementtinkerdescrips=bondpoltypeclassestoelementtinkerdescrips.copy()
        bondprms,bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips=RemoveOldParametersKeepNewParameters(poltype,bondprms,newbondindicestopoltypeclasses,'bond',bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips)
        strbndprms,anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips=RemoveOldParametersKeepNewParameters(poltype,strbndprms,newstrbndindicestopoltypeclasses,'strbnd',anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips)
        bondprms,bondpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,bondprms,bondindicestoextsmarts,bondsmartsatomordertoparameters,'bond',bondindicestoextsmartsatomorder)
        angleprms,anglepoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,angleprms,angleindicestoextsmarts,anglesmartsatomordertoparameters,'angle',angleindicestoextsmartsatomorder)
        strbndprms,strbndpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,strbndprms,angleindicestoextsmarts,strbndsmartsatomordertoparameters,'strbnd',angleindicestoextsmartsatomorder)
        torsionprms,torsionpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,torsionprms,torsionindicestoextsmarts,torsionsmartsatomordertoparameters,'torsion',torsionindicestoextsmartsatomorder)
        opbendprms,opbendpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,opbendprms,bondindicestoextsmarts,bondsmartsatomordertoparameters,'opbend',bondindicestoextsmartsatomorder)
        vdwprms,vdwpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,vdwprms,vdwindicestoextsmarts,vdwsmartsatomordertoparameters,'vdw',vdwindicestoextsmartsatomorder)
        tortorprms=[]
        tortorprms,tortorpoltypeclassestosmartsatomordersext=AddExternalDatabaseSMARTSMatchParameters(poltype,tortorprms,tortorindicestoextsmarts,tortorsmartsatomordertoparameters,'tortors',tortorindicestoextsmartsatomorders,tortorsmartsatomordertogrid)
        anglemissingtinkerclassestopoltypeclasses=RemovePoltypeClassesFromNewMatches(poltype,anglemissingtinkerclassestopoltypeclasses,anglepoltypeclassestoparametersmartsatomorders)
        bondmissingtinkerclassestopoltypeclasses=RemovePoltypeClassesFromNewMatches(poltype,bondmissingtinkerclassestopoltypeclasses,bondpoltypeclassestoparametersmartsatomorders)
        missinganglepoltypeclasses=ExtractMissingPoltypeClasses(poltype,anglemissingtinkerclassestopoltypeclasses)
        missingbondpoltypeclasses=ExtractMissingPoltypeClasses(poltype,bondmissingtinkerclassestopoltypeclasses)
        strbndprms=ZeroOutMissingStrbnd(poltype,anglemissingtinkerclassestopoltypeclasses,strbndprms)
        angleprms=AssignAngleGuessParameters(poltype,anglemissingtinkerclassestopoltypeclasses,angleprms,indextoneighbidxs)
        bondprms=AssignBondGuessParameters(poltype,bondmissingtinkerclassestopoltypeclasses,bondprms,indextoneighbidxs)
        angleprms.extend(newangleprms)
        bondprms.extend(newbondprms)
        strbndprms.extend(newstrbndprms)
        opbendprms.extend(newopbendprms)
        angleprms=ModifyAngleKeywords(poltype,angleprms,planarangletinkerclassestopoltypeclasses)
        bondlistbabel=ConvertToBabelList(poltype,listofbondsforprm)
        anglelistbabel=ConvertToBabelList(poltype,listofanglesforprm)
        bondprms=AddOptimizedBondLengths(poltype,optmol,bondprms,bondlistbabel)
        angleprms=AddOptimizedAngleLengths(poltype,optmol,angleprms,anglelistbabel)
        
        torsionsmissingtinkerclassestopoltypeclasses=MergeDictionaries(poltype,torsionszeroouttinkerclassestopoltypeclasses,torsionsmissingtinkerclassestopoltypeclasses)
        extratormissingtinkerclassestopoltypeclasses=GenerateTinkerClassesToPoltypeClasses(poltype,removedtorsionindices) 
        torsionsmissingtinkerclassestopoltypeclasses=MergeDictionaries(poltype,extratormissingtinkerclassestopoltypeclasses,torsionsmissingtinkerclassestopoltypeclasses)

        torsionprms=ZeroOutMissingTorsions(poltype,torsionsmissingtinkerclassestopoltypeclasses,torsionprms)
        torsionprms,arotorsionlinetodescrips=DefaultAromaticMissingTorsions(poltype,arotorsionsmissingtinkerclassestopoltypeclasses,partialarotorsionsmissingtinkerclassestopoltypeclasses,torsionprms,mol)
        torsionkeystringtoparameters=GrabTorsionParameterCoefficients(poltype,torsionprms)
        potentialmissingopbendprmtypes=FindPotentialMissingParameterTypes(poltype,opbendprms,planarbondtinkerclassestopoltypeclasses)
        potentialmissingopbendprmindices=ConvertPoltypeClassesToIndices(poltype,potentialmissingopbendprmtypes)
        potentialmissingopbendprmindices=FilterIndices(poltype,potentialmissingopbendprmindices,planarbonds)
        missingopbendprmindices=CheckIfParametersExist(poltype,potentialmissingopbendprmindices,opbendprms)
        torsionsmissing=ConvertToPoltypeClasses(poltype,torsionsmissing)
        missingvdwtypes=[poltype.idxtosymclass[i] for i in missingvdwatomindices]
        defaultvalues=None
        if len(missingopbendprmindices)!=0:
            newopbendprms,defaultvalues=DefaultOPBendParameters(poltype,missingopbendprmindices,mol,opbendbondindicestotrigonalcenterbools)
            opbendprms.extend(newopbendprms)
        opbendprms=AddTrivialOPBendForAmine(poltype,opbendprms,mol)
        if poltype.isfragjob==True and len(poltype.onlyrotbndslist)!=0:
            vdwprms,parentvdwtransferinfo=GrabVdwParametersFromParent(poltype,vdwprms)
        polarprmstotransferinfo=MapParameterLineToTransferInfo(poltype,newpolarprms,{},{},{},{},newpolarpoltypecommentstocomments,newpolarpoltypecommentstosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
        vdwprmstotransferinfo=MapParameterLineToTransferInfo(poltype,vdwprms,atompoltypeclasstoparametersmartsatomorder,atompoltypeclasstosmartsatomorder,atompoltypeclassestoelementtinkerdescrip,vdwpoltypeclassestosmartsatomordersext,{},{},arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
        if poltype.forcefield=='AMOEBA+':
            amoebaplusvdwprmstotransferinfo=MapParameterLineToTransferInfo(poltype,newvdwprms,{},{},{},{},newvdwpoltypecommentstocomments,newvdwpoltypecommentstosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
            ctprmstotransferinfo=MapParameterLineToTransferInfo(poltype,ctprms,{},{},{},{},ctpoltypecommentstocomments,ctpoltypecommentstosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
            cpprmstotransferinfo=MapParameterLineToTransferInfo(poltype,cpprms,{},{},{},{},cppoltypecommentstocomments,cppoltypecommentstosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)

            bondcfprmstotransferinfo=MapParameterLineToTransferInfo(poltype,bondcfprms,{},{},{},{},bondcfpoltypecommentstocomments,bondcfpoltypecommentstosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
            anglecfprmstotransferinfo=MapParameterLineToTransferInfo(poltype,anglecfprms,{},{},{},{},anglecfpoltypecommentstocomments,anglecfpoltypecommentstosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)


        else:
            amoebaplusvdwprmstotransferinfo={}
            ctprmstotransferinfo={}
            cpprmstotransferinfo={}
            bondcfprmstotransferinfo={}
            anglecfprmstotransferinfo={}
        tortorprmstotransferinfo=MapParameterLineToTransferInfo(poltype,tortorprms,{},{},{},tortorpoltypeclassestosmartsatomordersext,{},{},arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses,defaultvalues=None,keyword='tortors')
        bondprmstotransferinfo=MapParameterLineToTransferInfo(poltype,bondprms,bondpoltypeclassestoparametersmartsatomorders,bondpoltypeclassestosmartsatomorders,bondpoltypeclassestoelementtinkerdescrips,bondpoltypeclassestosmartsatomordersext,newbondpoltypeclassestocomments,newbondpoltypeclassestosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
        opbendprmstotransferinfo=MapParameterLineToTransferInfo(poltype,opbendprms,blankbondpoltypeclassestoparametersmartsatomorders,opbendpoltypeclassestosmartsatomorders,opbendpoltypeclassestoelementtinkerdescrips,opbendpoltypeclassestosmartsatomordersext,newopbendpoltypeclassestocomments,newopbendpoltypeclassestosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses,defaultvalues)
        
        angleprmstotransferinfo=MapParameterLineToTransferInfo(poltype,angleprms,anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips,anglepoltypeclassestosmartsatomordersext,newanglepoltypeclassestocomments,newanglepoltypeclassestosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
        strbndprmstotransferinfo=MapParameterLineToTransferInfo(poltype,strbndprms,anglepoltypeclassestoparametersmartsatomorders,anglepoltypeclassestosmartsatomorders,anglepoltypeclassestoelementtinkerdescrips,strbndpoltypeclassestosmartsatomordersext,newstrbndpoltypeclassestocomments,newstrbndpoltypeclassestosmartslist,arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
        torsionprmstotransferinfo=MapParameterLineToTransferInfo(poltype,torsionprms,torsionpoltypeclassestoparametersmartsatomorders,torsionpoltypeclassestosmartsatomorders,torsionpoltypeclassestoelementtinkerdescrips,torsionpoltypeclassestosmartsatomordersext,{},{},arotorsionlinetodescrips,missingvdwtypes,torsionsmissing,missingbondpoltypeclasses,missinganglepoltypeclasses)
        if poltype.isfragjob==True and len(poltype.onlyrotbndslist)!=0:
            vdwprmstotransferinfo=AddParentVdwTransferInfo(poltype,parentvdwtransferinfo,vdwprmstotransferinfo)
        WriteOutList(poltype,torsionsmissing,poltype.torsionsmissingfilename)
        WriteDictionaryToFile(poltype,torsionkeystringtoparameters,poltype.torsionprmguessfilename)
        WriteOutList(poltype,missingvdwatomindices,poltype.vdwmissingfilename)
        WriteOutList(poltype,tortorsmissing,poltype.tortormissingfilename)
        return bondprmstotransferinfo,angleprmstotransferinfo,torsionprmstotransferinfo,strbndprmstotransferinfo,opbendprmstotransferinfo,vdwprmstotransferinfo,polarprmstotransferinfo,torsionsmissing,torsionkeystringtoparameters,missingvdwatomindices,soluteprms,amoebaplusvdwprmstotransferinfo,ctprmstotransferinfo,cpprmstotransferinfo,bondcfprmstotransferinfo,anglecfprmstotransferinfo,tortorprmstotransferinfo,tortorsmissing
