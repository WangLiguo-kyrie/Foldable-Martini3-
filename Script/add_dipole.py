#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 13:31:48 2022

@author: Wang chris
used for add dipolement on protein Backbone as in ProMPT
"""

import vermouth
from vermouth.forcefield import ForceField
import pandas as pd
import argparse

if __name__ == '__main__':
    
    ff = ForceField('martini3001')

    parser = argparse.ArgumentParser()
    parser.add_argument("molecule", help="name of peptide to process", type = str)
    args = parser.parse_args()
    
    file = '%s.itp' %args.molecule
    
    l = []
    with open(file, 'r') as f:
        for line in f.readlines():
            l.append(line)
            
    mol = vermouth.gmx.read_itp(l, ff)
    
    molname = args.molecule
    
    df = pd.DataFrame(dict(ff.blocks[molname].nodes))
    
    BB_inds = df.loc['index'][df.loc['atomname'] == 'BB'].values
    resids = df.loc['resid'][df.loc['atomname'] == 'BB'].values
    resnames= df.loc['resname'][df.loc['atomname'] == 'BB'].values
    
    print(df)
    num_res=len(BB_inds)
    print('Reside number:{}'.format(num_res))
    num_atom=df.loc['index'].values[-1]
    print('Atom number:{}'.format(num_atom))
    print(resnames)
    
    # add dummy bead dipole in [atoms]
    # Specify Pro AND Gly. Dipole charge and LJ protential were both removed for Pro ((LJ for negative potetnial remained)),
    # only LJ potential was removed for Gly to reduce these two residue's ss formation propensity. 
    for i in range(1, 1+num_res):
        print(i)
        resname_i=resnames[i-1]
        print(resname_i)
        if (resname_i == 'PRO'):
            ff.blocks[molname].add_node(i+num_atom-1,
                                    index=i+num_atom,
                                    atomname='UP',
                                    atype='UPP',
                                    resname=resnames[i-1],
                                    resid=resids[i-1],
                                    charge_group=i,
                                    charge=0.0,
                                    mass = 24)
        elif (resname_i == 'GLY'):
            ff.blocks[molname].add_node(i+num_atom-1,
                                    index=i+num_atom,
                                    atomname='UP',
                                    atype='UPG',
                                    resname=resnames[i-1],
                                    resid=resids[i-1],
                                    charge_group=i,
                                    charge=0.34,
                                    mass = 24)
        elif ((resname_i == 'ARG')|(resname_i == 'ALA')|(resname_i == 'LEU')|(resname_i == 'LYS')|(resname_i == 'MET')):   
            ff.blocks[molname].add_node(i+num_atom-1,
                                    index=i+num_atom,
                                    atomname='UP',
                                    atype='UP1',
                                    resname=resnames[i-1],
                                    resid=resids[i-1],
                                    charge_group=i,
                                    charge=0.34,
                                    mass = 24)                                         
        elif ((resname_i == 'GLU')|(resname_i == 'GLN')|(resname_i == 'ILE')|(resname_i == 'PHE')|(resname_i == 'SER')|(resname_i == 'TYR')|(resname_i == 'TRP')):   
            ff.blocks[molname].add_node(i+num_atom-1,
                                    index=i+num_atom,
                                    atomname='UP',
                                    atype='UP2',
                                    resname=resnames[i-1],
                                    resid=resids[i-1],
                                    charge_group=i,
                                    charge=0.34,
                                    mass = 24) 
        else:   
            ff.blocks[molname].add_node(i+num_atom-1,
                                    index=i+num_atom,
                                    atomname='UP',
                                    atype='UP3',
                                    resname=resnames[i-1],
                                    resid=resids[i-1],
                                    charge_group=i,
                                    charge=0.34,
                                    mass = 24) 
                                    
    for i in range(1, 1+num_res):
        print(i)
        resname_i=resnames[i-1]
        print(resname_i)
        if (resname_i == 'PRO'):
            ff.blocks[molname].add_node(i+num_atom+num_res-1,
                                    index=i+num_atom+num_res,
                                    atomname='UN',
                                    atype='UN2',
                                    resname=resnames[i-1],
                                    resid=resids[i-1],
                                    charge_group=i,
                                    charge=0.0,
                                    mass = 24)
        elif (resname_i == 'GLY'):
            ff.blocks[molname].add_node(i+num_atom+num_res-1,
                                    index=i+num_atom+num_res,
                                    atomname='UN',
                                    atype='UNG',
                                    resname=resnames[i-1],
                                    resid=resids[i-1],
                                    charge_group=i,
                                    charge=-0.34,
                                    mass = 24)
        elif ((resname_i == 'ARG')|(resname_i == 'ALA')|(resname_i == 'LEU')|(resname_i == 'LYS')|(resname_i == 'MET')):     
            ff.blocks[molname].add_node(i+num_atom+num_res-1,
                                    index=i+num_atom+num_res,
                                    atomname='UN',
                                    atype='UN1',
                                    resname=resnames[i-1],
                                    resid=resids[i-1],
                                    charge_group=i,
                                    charge=-0.34,
                                    mass = 24)
        elif ((resname_i == 'GLU')|(resname_i == 'GLN')|(resname_i == 'ILE')|(resname_i == 'PHE')|(resname_i == 'SER')|(resname_i == 'TYR')|(resname_i == 'TRP')):      
            ff.blocks[molname].add_node(i+num_atom+num_res-1,
                                    index=i+num_atom+num_res,
                                    atomname='UN',
                                    atype='UN2',
                                    resname=resnames[i-1],
                                    resid=resids[i-1],
                                    charge_group=i,
                                    charge=-0.34,
                                    mass = 24)                                    
        else:      
            ff.blocks[molname].add_node(i+num_atom+num_res-1,
                                    index=i+num_atom+num_res,
                                    atomname='UN',
                                    atype='UN3',
                                    resname=resnames[i-1],
                                    resid=resids[i-1],
                                    charge_group=i,
                                    charge=-0.34,
                                    mass = 24)   
                                                                        
    df = pd.DataFrame(dict(ff.blocks[molname].nodes))
    print(df)
    # add backbone-dummy bead dipole bond in [bonds]
    for i in range(1, 1+num_res):
        print(i)
        ff.blocks[molname].add_interaction('bonds',
                                       atoms = [BB_inds[i-1]-1, i+num_atom-1],
                                       parameters = ['1', '0.140', '10000'],
                                       meta={'comment':'BB-UP stiff bond'}
                                       )
        ff.blocks[molname].add_interaction('bonds',
                                       atoms = [BB_inds[i-1]-1, i+num_atom+num_res-1],
                                       parameters = ['1', '0.140', '10000'],
                                       meta={'comment':'BB-UN stiff bond'}
                                       )    
                                       
    # add backbone-dummy bead dipole bond constraint in [constraints]
    for i in range(1, 1+num_res):
        print(i)
        ff.blocks[molname].add_interaction('constraints',
                                       atoms = [BB_inds[i-1]-1, i+num_atom-1],
                                       parameters = ['1', '0.140'],
                                       meta={'comment':'BB-UP constraint'}
                                       )
        ff.blocks[molname].add_interaction('constraints',
                                       atoms = [BB_inds[i-1]-1, i+num_atom+num_res-1],
                                       parameters = ['1', '0.140'],
                                       meta={'comment':'BB-UN constraint'}
                                       )    
                                       
    # add UN/UP-BB-BB and BB-BB-UN/UP angles in [angles] like BBS/SBB, keeping dipole perpendicular to BB-BB peptide bond
    for i in range(1, num_res):
        print(i)                                       
        ff.blocks[molname].add_interaction('angles',
                                       atoms = [i+num_atom-1,BB_inds[i-1]-1,BB_inds[i]-1],
                                       parameters = ['10', '100', '15'],
                                       meta={'comment':'UP-BB-BB'}
                                       )                                       
        ff.blocks[molname].add_interaction('angles',
                                       atoms = [i+num_atom+num_res-1,BB_inds[i-1]-1,BB_inds[i]-1],
                                       parameters = ['10', '100', '15'],
                                       meta={'comment':'UN-BB-BB'}
                                       )                                         
    for i in range(2, num_res+1):
        print(i)                                       
        ff.blocks[molname].add_interaction('angles',
                                       atoms = [BB_inds[i-2]-1,BB_inds[i-1]-1,i+num_atom-1],
                                       parameters = ['10', '100', '15'],
                                       meta={'comment':'BB-BB-UP'}
                                       )                                       
        ff.blocks[molname].add_interaction('angles',
                                       atoms = [BB_inds[i-2]-1,BB_inds[i-1]-1,i+num_atom+num_res-1],
                                       parameters = ['10', '100', '15'],
                                       meta={'comment':'BB-BB-UN'}
                                       )  
                                       
    # define orientation restrain occurs in the peptide bond planar, UN in resid-i and UP in next resid-i+1 to be in opposite direction
    # NO restrain for the orientation of UP/UN in the same residue
    for i in range(1, num_res):
        print(i)     
        ff.blocks[molname].add_interaction('dihedrals',
                                       atoms = [i+num_atom+num_res-1,BB_inds[i-1]-1,BB_inds[i]-1,i+num_atom+1-1],
                                       parameters = ['1', '0', '100', '1'],
                                       meta={'comment':'UN-BB-BB-UP planar'}
                                       )                                                                                                                                                                                                                

    # add dipole-dipole Exclusions    
    # caution: the atom_index in parameter option doesn't need -1
    for i in range(1, 1+num_res):
        print(i)
        ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [i+num_atom-1],
                                       parameters = [i+num_atom+num_res],
                                       meta={'comment':'UP-UN in same residue'}
                                       )
    for i in range(1, num_res):
        print(i)
        ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [i+num_atom-1],parameters = [i+num_atom+1],
                                       meta={'comment':'UP-UP in 1-2 residue'}
                                       )
        ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [i+num_atom-1],parameters =[i+num_atom+num_res+1],
                                       meta={'comment':'UP-UN in 1-2 residue'}
                                       )
        ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [i+num_atom+num_res-1],parameters =[ i+num_atom+1],
                                       meta={'comment':'UN-UP in 1-2 residue'}
                                       )
        ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [i+num_atom+num_res-1],parameters =[ i+num_atom+num_res+1],
                                       meta={'comment':'UN-UN in 1-2 residue'}
                                       )
    # add dipole-charged sidechain Exclusions  
    for i in range(1, 1+num_res):
        print(i)
        resname_i=resnames[i-1]
        print(resname_i)
        if ((resname_i == 'ARG')|(resname_i == 'LYS')):
            ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [i+num_atom-1],parameters =[ BB_inds[i-1]+2],
                                       meta={'comment':'UP-charged sidechain in same residue'}
                                       )    
            ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [i+num_atom+num_res-1],parameters =[ BB_inds[i-1]+2],
                                       meta={'comment':'UN-charged sidechain in same residue'}
                                       )    
        elif ((resname_i == 'GLU')|(resname_i == 'ASP')):
            ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [i+num_atom-1],parameters =[ BB_inds[i-1]+1],
                                       meta={'comment':'UP-charged sidechain in same residue'}
                                       )    
            ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [i+num_atom+num_res-1],parameters =[ BB_inds[i-1]+1],
                                       meta={'comment':'UN-charged sidechain in same residue'}
                                       )                                          
    # add charged terminal and backbone diploe  Exclusions
    ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [0],parameters =[ num_atom+1],
                                       meta={'comment':'Nter-UP in same residue'}
                                       )
    ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [0],parameters =[ num_atom+num_res+1],
                                       meta={'comment':'Nter-UN in same residue'}
                                       )
    ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [0],parameters =[ num_atom+2],
                                       meta={'comment':'Nter-UP in 1-2 residue'}
                                       )
    ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [0],parameters =[ num_atom+num_res+2],
                                       meta={'comment':'Nter-UN in 1-2 residue'}
                                       )  
                                       
    ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [BB_inds[num_res-1]-1],parameters =[ num_atom+num_res],
                                       meta={'comment':'Cter-UP in same residue'}
                                       )
    ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [BB_inds[num_res-1]-1],parameters =[ num_atom+num_res*2],
                                       meta={'comment':'Cter-UN in same residue'}
                                       )  
    ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [BB_inds[num_res-1]-1],parameters =[ num_atom+num_res-1],
                                       meta={'comment':'Cter-UP in 1-2 residue'}
                                       )
    ff.blocks[molname].add_interaction('exclusions',
                                       atoms = [BB_inds[num_res-1]-1],parameters =[ num_atom+num_res*2-1],
                                       meta={'comment':'Cter-UN in 1-2 residue'}
                                       )                                       
                                                                                                                 

    ff.blocks[molname].meta['moltype'] = molname
    
    mol_out = ff.blocks[molname].to_molecule()    
    mol_out.meta['moltype'] = molname
    
    header_lines = ['itp modified by Wang Chris used for add dipolement on protein Backbone']
    
    with open('mol.itp', 'w') as outfile:
        vermouth.gmx.write_molecule_itp(mol_out, outfile=outfile, header=header_lines)
    
