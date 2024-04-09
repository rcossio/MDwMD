#!/bin/bash

# Useful
# MDP options: https://manual.gromacs.org/2023.2/user-guide/mdp-options.html
# File formats: https://manual.gromacs.org/2023.2/reference-manual/file-formats.html

# Details
# Force field setup: https://manual.gromacs.org/2023.2/user-guide/force-fields.html#charmm

# Convert to gromacs format
gmx pdb2gmx -f my_proteina.pdb -o protein_processed.gro -p topol.top -i posre.itp -posrefc 1000 -ff charmm36-jul2022 -water tip3p -ignh -nobackup
cp posre.itp posre.itp.bak


# Genbox and solvate
gmx editconf -f protein_processed.gro -o protein_newbox.gro -c -d 1.5 -bt cubic -nobackup
gmx solvate -cp protein_newbox.gro -cs spc216.gro -o protein_solv.gro -p topol.top -nobackup


# Add ions
touch ions.mdp
gmx grompp -f ions.mdp -c protein_solv.gro -o ions.tpr -p topol.top
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -conc 0.15 -neutral -nobackup <<< "SOL"


# Minimization
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -nb gpu
gmx energy -f em.edr -o em_potential.xvg -xvg none -nobackup <<< $'Potential\n0\n'


# NVT restrained
gmx grompp -f nvt_r.mdp -c em.gro -r em.gro -p topol.top -o nvt_r.tpr
gmx mdrun -deffnm nvt_r -nb gpu
gmx energy -f nvt_r.edr -o nvt_r_temperature.xvg -xvg none -nobackup <<< $'Temperature\n0\n'


# Reducing restraints
gmx grompp -f npt_r1000.mdp -c nvt_r.gro -t nvt_r.cpt -r em.gro -p topol.top -o npt_r1000.tpr
gmx mdrun -deffnm npt_r1000 -nb gpu
gmx energy -f npt_r1000.edr -o npt_r1000.xvg -xvg none -nobackup <<< $'Pressure\nDensity\nPosition-Rest\nPotential\n0\n'

sed -e "s/1000/200/g" posre.itp.bak > posre.itp
gmx grompp -f npt_r200.mdp -c npt_r1000.gro -t npt_r1000.cpt -r em.gro -p topol.top -o npt_r200.tpr
gmx mdrun -deffnm npt_r200 -nb gpu
gmx energy -f npt_r200.edr -o npt_r200.xvg -xvg none -nobackup <<< $'Pressure\nDensity\nPosition-Rest\nPotential\n0\n'

sed -e "s/1000/50/g" posre.itp.bak > posre.itp
gmx grompp -f npt_r50.mdp -c npt_r200.gro -t npt_r200.cpt -r em.gro -p topol.top -o npt_r50.tpr
gmx mdrun -deffnm npt_r50 -nb gpu
gmx energy -f npt_r50.edr -o npt_r50.xvg -xvg none -nobackup <<< $'Pressure\nDensity\nPosition-Rest\nPotential\n0\n'

sed -e "s/1000/10/g" posre.itp.bak > posre.itp
gmx grompp -f npt_r10.mdp -c npt_r50.gro -t npt_r50.cpt -r em.gro -p topol.top -o npt_r10.tpr
gmx mdrun -deffnm npt_r10 -nb gpu
gmx energy -f npt_r10.edr -o npt_r10.xvg -xvg none -nobackup <<< $'Pressure\nDensity\nPosition-Rest\nPotential\n0\n'

sed -e "s/1000/2/g" posre.itp.bak > posre.itp
gmx grompp -f npt_r2.mdp -c npt_r10.gro -t npt_r10.cpt -r em.gro -p topol.top -o npt_r2.tpr
gmx mdrun -deffnm npt_r2 -nb gpu
gmx energy -f npt_r2.edr -o npt_r2.xvg -xvg none -nobackup <<< $'Pressure\nDensity\nPosition-Rest\nPotential\n0\n'


# Free NPT MD
gmx grompp -f npt_f.mdp -c npt_r2.gro -t npt_r2.cpt -p topol.top -o npt_f.tpr
gmx mdrun -deffnm npt_f -nb gpu
gmx energy -f npt_f.edr -o npt_f.xvg -xvg none -nobackup <<< $'Pressure\nDensity\n0\n'


# Free NVT MD
gmx grompp -f nvt_f.mdp -c npt_f.gro -t npt_f.cpt -p topol.top -o nvt_f.tpr 
gmx mdrun -deffnm nvt_f -nb gpu
gmx energy -f nvt_f.edr -o nvt_f.xvg -xvg none -nobackup <<< $'Pressure\n0\n'


