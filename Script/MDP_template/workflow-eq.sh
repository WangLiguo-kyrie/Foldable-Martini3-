#!/bin/bash
gmx grompp -f martini_em.mdp -c solv-ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -v -nt 10 -tableb table_d1.xvg
gmx grompp -f martini_nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -v -nt 10 -tableb table_d1.xvg
gmx grompp -f martini_npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -v -nt 10 -tableb table_d1.xvg

gmx grompp -f martini_npt2.mdp -c npt.gro -p topol.top -o npt2.tpr
gmx mdrun -deffnm npt2 -v -nt 10 -tableb table_d1.xvg
gmx grompp -f martini_npt3.mdp -c npt2.gro -p topol.top -o npt3.tpr
gmx mdrun -deffnm npt3 -v -nt 10 -tableb table_d1.xvg
gmx grompp -f martini_npt4.mdp -c npt3.gro -p topol.top -o npt4.tpr
gmx mdrun -deffnm npt4 -v -nt 10 -tableb table_d1.xvg
gmx grompp -f martini_npt5.mdp -c npt4.gro -p topol.top -o npt5.tpr
gmx mdrun -deffnm npt5 -v -nt 10 -tableb table_d1.xvg
gmx grompp -f martini_run.mdp -c npt5.gro -p topol.top -o production.tpr -maxwarn 1
gmx mdrun -deffnm production -v -nt 10 -tableb table_d1.xvg
