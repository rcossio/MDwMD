#!/bin/bash

# Useful
# MDP options: https://manual.gromacs.org/2023.2/user-guide/mdp-options.html
# File formats: https://manual.gromacs.org/2023.2/reference-manual/file-formats.html

# Details
# Force field setup: https://manual.gromacs.org/2023.2/user-guide/force-fields.html#charmm

# Convert to gromacs format
gmx pdb2gmx -f my_proteina.pdb -o protein_processed.gro -p topol.top -i posre.itp -posrefc 1000 -ff charmm36-jul2022 -water tip3p -ignh -nobackup


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


# NVT
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -nb gpu
gmx energy -f nvt.edr -o nvt_temperature.xvg -xvg none -nobackup <<< $'Temperature\n0\n'


# NPT
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r em.gro -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -nb gpu
gmx energy -f npt.edr -o npt.xvg -xvg none -nobackup <<< $'Pressure\nDensity\nPosition-Rest\nPotential\n0\n'


# Reducing restraints
sed -i 's/1000/200/g' posre.itp 
gmx grompp -f rst200.mdp -c npt.gro -t npt.cpt -r em.gro -p topol.top -o rst_200.tpr
gmx mdrun -deffnm rst_200 -nb gpu
gmx energy -f rst_200.edr -o rst_200.xvg -xvg none -nobackup <<< $'Pressure\nDensity\nPosition-Rest\nPotential\n0\n'

sed -i 's/200/50/g' posre.itp 
gmx grompp -f rst50.mdp -c rst_200.gro -t rst_200.cpt -r em.gro -p topol.top -o rst_50.tpr
gmx mdrun -deffnm rst_50 -nb gpu
gmx energy -f rst_50.edr -o rst_50.xvg -xvg none -nobackup <<< $'Pressure\nDensity\nPosition-Rest\nPotential\n0\n'

sed -i 's/50/10/g' posre.itp 
gmx grompp -f rst10.mdp -c rst_50.gro -t rst_50.cpt -r em.gro -p topol.top -o rst_10.tpr
gmx mdrun -deffnm rst_10 -nb gpu
gmx energy -f rst_10.edr -o rst_10.xvg -xvg none -nobackup <<< $'Pressure\nDensity\nPosition-Rest\nPotential\n0\n'

gmx grompp -f rst0.mdp -c rst_10.gro -t rst_10.cpt -p topol.top -o rst_0.tpr -maxwarn 1 #Warning with center of mass removal
gmx mdrun -deffnm rst_0 -nb gpu
gmx energy -f rst_0.edr -o rst_0.xvg -xvg none -nobackup <<< $'Pressure\nDensity\n0\n'

# Free NVT MD
gmx grompp -f md.mdp -c rst_0.gro -t rst_0.cpt -p topol.top -o md.tpr  -maxwarn 1 #Warning with center of mass removal
gmx mdrun -deffnm md -nb gpu
gmx energy -f md.edr -o md.xvg -xvg none -nobackup <<< $'Pressure\n0\n'


