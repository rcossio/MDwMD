#!/bin/bash

# Useful
# MDP options: https://manual.gromacs.org/2023.2/user-guide/mdp-options.html
# File formats: https://manual.gromacs.org/2023.2/reference-manual/file-formats.html

# Details
# Force field setup: https://manual.gromacs.org/2023.2/user-guide/force-fields.html#charmm


# Pasar a formato gromacs
gmx pdb2gmx -f tu_proteina.pdb -o protein_processed.gro -p topol.top -i posre.itp -posrefc 1000 -ff charmm36-jul2022 -water tip3p -ignh -nobackup


# Crear caja de simulación
gmx editconf -f protein_processed.gro -o protein_newbox.gro -c -d 1.0 -bt cubic -nobackup

# Solvatar el sistema
gmx solvate -cp protein_newbox.gro -cs spc216.gro -o protein_solv.gro -p topol.top -nobackup

# Agregar iones
touch ions.mdp

gmx grompp -f ions.mdp -c protein_solv.gro -o ions.tpr -p topol.top

gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -conc 0.15 -neutral -nobackup # Poner 13

# Minimización
cat > em.mdp << EOF
integrator      = steep     ; Steepest descent minimization
emtol           = 100.0     ; Tolerance for energy minimization
emstep          = 0.01      ; Minimization step size
nsteps          = 10000     ; Max steps for energy minimization
nstlist         = 10        ; Update neighbor list every step
rlist           = 1.2       ; Cut-off for making neighbor list (short range forces)
cutoff-scheme   = Verlet    ; Use Verlet list for neighbor search
ns-type         = grid      ; Grid search for neighbor list
coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb        = 1.2       ; Cut-off for Coulomb interactions in nm
rvdw            = 1.2       ; Cut-off for Van der Waals interactions in nm
constraints     = h-bonds 
constraint-algorithm = LINCS 
DispCorr        = no
EOF

gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr

gmx mdrun -v -deffnm em -nb gpu

gmx energy -f em.edr -o em_potential.xvg -xvg none -nobackup # Poner 11 0

# Equilibración de temperatura
cat > nvt.mdp << EOF
; Apply position restraints
define                  = -DPOSRES

; Run parameters
integrator              = md 
dt                      = 0.002     
nsteps                  = 1000000

; Output control
nstxout-compressed      = 10000   ; coordinates
nstvout                 = 10000   ; velocities
nstcalcenergy 		    = 100    ; energie calculation
nstenergy               = 10000   ; energy saving
nstlog                  = 10000   ; log file

; Non-bonded parameters
cutoff-scheme           = Verlet 
nstlist                 = 20  
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2

; Electrostatics
coulombtype             = PME       
rcoulomb                = 1.2      
DispCorr                = no

; Temperature coupling
tcoupl                  = V-rescale                 
tc-grps                 = Protein Non-Protein
tau-t                   = 0.1 0.1                   
ref-t                   = 293 293

; Bond constraints
constraints             = h-bonds
constraint-algorithm    = LINCS 
continuation            = no 

; Velocity generation
gen-vel                 = yes       
gen-temp                = 293    
gen-seed                = -1

; Center of mass
comm-mode               = None
EOF

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

gmx mdrun -v -deffnm nvt -nb gpu

gmx energy -f nvt.edr -o nvt_temperature.xvg -xvg none -nobackup # Poner 16 0


# Equilibración de presión
cat > npt.mdp <<EOF
; Apply position restraints
define                  = -DPOSRES

; Run parameters
integrator              = md      
dt                      = 0.002
nsteps                  = 1000000

; Output control
nstxout-compressed      = 10000
nstvout                 = 10000
nstfout	          	    = 10000
nstcalcenergy 		    = 100
nstenergy               = 10000     
nstlog                  = 10000       

; Non-bonded parameters
cutoff-scheme           = Verlet
nstlist                 = 20      
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2      

; Electrostatics
coulombtype             = PME    
rcoulomb                = 1.2
DispCorr                = no

; Temperature coupling
tcoupl                  = V-rescale               
tc-grps                 = System
tau-t                   = 0.1
ref-t                   = 293

; Pressure coupling
pcoupl                  = C-rescale                   
pcoupltype              = isotropic                 
tau-p                   = 2.0                          
ref-p                   = 1.0                       
compressibility         = 4.5e-5                     
refcoord-scaling        = all

; Bond constraints
continuation            = yes       
constraints             = h-bonds 
constraint-algorithm    = LINCS

; Velocity generation
gen-vel                 = no   

; Center of mass
comm-mode               = None
EOF

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r em.gro -p topol.top -o npt.tpr

gmx mdrun -v -deffnm npt -nb gpu

gmx energy -f npt.edr -o npt_pressure.xvg -xvg none -nobackup # Poner 17 0

gmx energy -f npt.edr -o npt_density.xvg -xvg none -nobackup # Poner 23 0

gmx energy -f npt.edr -o rst_1000.xvg -xvg none -nobackup # Poner 11 12 0


# Eliminación de restricciones
cp npt.mdp rst.mdp

sed -i 's/1000/200/g' posre.itp 

gmx grompp -f rst.mdp -c npt.gro -t npt.cpt -r em.gro -p topol.top -o rst_200.tpr

gmx mdrun -v -deffnm rst_200 -nb gpu

gmx energy -f rst_200.edr -o rst_200.xvg -xvg none -nobackup # Poner 11 12 0


sed -i 's/200/50/g' posre.itp 
gmx grompp -f rst.mdp -c rst_200.gro -t rst_200.cpt -r em.gro -p topol.top -o rst_50.tpr

gmx mdrun -v -deffnm rst_50 -nb gpu

gmx energy -f rst_50.edr -o rst_50.xvg -xvg none -nobackup # Poner 11 12 0


sed -i 's/50/10/g' posre.itp 
gmx grompp -f rst.mdp -c rst_50.gro -t rst_50.cpt -r em.gro -p topol.top -o rst_10.tpr

gmx mdrun -v -deffnm rst_10 -nb gpu

gmx energy -f rst_10.edr -o rst_10.xvg -xvg none -nobackup # Poner 11 12 0