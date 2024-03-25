#!/bin/bash

# Useful
# MDP options: https://manual.gromacs.org/2023.2/user-guide/mdp-options.html
# File formats: https://manual.gromacs.org/2023.2/reference-manual/file-formats.html

# Details
# Force field setup: https://manual.gromacs.org/2023.2/user-guide/force-fields.html#charmm


export PATH="/opt/gromacs/2023.2/bin/:$PATH"

input="5PTI"

# 1. Convert PDB to GROMACS format, selecting force field and water model
gmx pdb2gmx -f "5PTI.pdb" -o protein_processed.gro -p topol.top -ff charmm36-jul2022 -water tip3p -ignh -nobackup
# + files: posre.itp  protein_processed.gro  topol.top

# 2. Define the simulation box
gmx editconf -f protein_processed.gro -o protein_newbox.gro -c -d 1.0 -bt cubic -nobackup
# + files: protein_newbox.gro
# edits: topol.top 

# 3. Solvate the system
gmx solvate -cp protein_newbox.gro -cs spc216.gro -o protein_solv.gro -p topol.top -nobackup
# + files: protein_solv.gro
# edits: topol.top

# 4. Add ions to neutralize the system 
cat > ions.mdp << EOF
integrator      = steep     ; Algorithm (steepest descent minimization)
emtol           = 1000.0    ; Stop minimization upon maximum force
emstep          = 0.01      ; Energy step size
nsteps          = 50000     ; Maximum number of steps to perform
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme	 = Verlet    ; Buffered neighbor searching 
rlist		       = 1.0		 ; Cut-off for making neighbor list (short range forces)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
EOF
# + files: ions.mdp

gmx grompp -f ions.mdp -c protein_solv.gro -o ions.tpr -p topol.top
# + files: mdout.mdp ions.tpr

gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -conc 0.15 -neutral -nobackup << EOF
13
EOF
# + files: solv_ions.gro
# edits: topol.top

# Generate group to restrain the protein backbone
gmx make_ndx -f solv_ions.gro -o index.ndx << EOF
q
EOF
# + files: index.ndx

# 6. Energy Minimization
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
constraints          = h-bonds 
constraint-algorithm = LINCS 
EOF 
# + files: em.mdp

gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
# + files: em.tpr mdout.mdp

gmx mdrun -v -deffnm em -nb gpu
# + files: em.gro em.trr em.log em.edr

#To check
gmx energy -f em.edr -o potential.xvg << EOF
11
0
EOF
# + files: potential.xvg

# 7. NVT Equilibration (with position restraints)
#tc-grps = Protein Non-Protein
cat > nvt.mdp << EOF
integrator              = md 
dt                      = 0.002  ; 2 fs timestep     
nsteps                  = 1000000; 2 ns simulation
nstxout-compressed      = 5000   ; save coordinates every 10 ps
nstvout                 = 5000   ; save velocities every 10 ps

nstcalcenergy 		      = 100    ; calculate energies every 0.2 ps
nstenergy               = 1000   ; save energies every 2 ps
nstlog                  = 1000   ; update log file every 2 ps

cutoff-scheme           = Verlet 
nstlist                 = 20  
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME       
rcoulomb                = 1.2      

tcoupl                  = V-rescale                 
tc-grps                 = System
tau-t                   = 0.1                   
ref-t                   = 300
        
continuation            = no 
constraints             = h-bonds
constraint-algorithm    = LINCS 

gen-vel                 = yes       
gen-temp                = 300       
gen-seed                = -1

define                  = -DPOSRES
EOF
# + files: nvt.mdp

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
# + files: nvt.tpr mdout.mdp

gmx mdrun -v -deffnm nvt -nb gpu
# + files: nvt.log nvt.trr nvt.edr nvt.xtc nvt.gro nvt.cpt

#To check
gmx energy -f nvt.edr -o temperature.xvg << EOF
16
0
EOF
# + files: temperature.xvg

# 8. NPT Equilibration (also with position restrains)
cat > npt.mdp <<EOF
integrator              = md      
dt                      = 0.002
nsteps                  = 1000000
nstxout-compressed      = 5000
nstvout                 = 5000
nstfout	          	   = 5000
nstcalcenergy 		      = 100
nstenergy               = 1000     
nstlog                  = 1000       


cutoff-scheme           = Verlet
nstlist                 = 20      
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       
coulombtype             = PME    
rcoulomb                = 1.2

tcoupl                  = V-rescale               
tc-grps                 = System 
tau_t                   = 0.1
ref_t                   = 300

pcoupl                  = C-rescale                   
pcoupltype              = isotropic                 
tau_p                   = 2.0                          
ref_p                   = 1.0                       
compressibility         = 4.5e-5                     
refcoord_scaling        = com

continuation            = yes       
constraints             = h-bonds 
constraint_algorithm    = LINCS

gen_vel                 = no   

define                  = -DPOSRES
EOF
# + files: npt.mdp

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
# + files: npt.tpr mdout.mdp

gmx mdrun -v -deffnm npt -nb gpu
# + files: npt.log npt.trr npt.xtc npt.edr npt.gro npt.cpt npt.tpr

#To check
gmx energy -f npt.edr -o pressure.xvg << EOF
17
0
EOF
# + files: pressure.xvg

gmx energy -f npt.edr -o density.xvg << EOF
23
0
EOF
# + files: density.xvg

# 9. Production MD
cat > md.mdp <<EOF
; Protein-ligand complex NVT equilibration 
define                  = -DPOSRES
integrator              = md 
dt                      = 0.002
nsteps                  = 25000000
nstxout-compressed      = 5000
nstvout                 = 5000
nstfout	          	   = 5000
nstcalcenergy 		      = 100
nstenergy               = 1000 
nstlog                  = 1000  
;
cutoff-scheme           = Verlet
; ns_type                 = grid
nstlist                 = 20  
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME       
rcoulomb                = 1.2      
; pme_order               = 4     
; fourierspacing          = 0.16  
;
tcoupl                  = V-rescale                 
tc-grps                 = System
tau_t                   = 0.1   0.1                   
ref_t                   = 300   300
;         
continuation            = no 
constraints             = h-bonds
constraint_algorithm    = LINCS 
; lincs_iter              = 1      
; lincs_order             = 4 
;
; pcoupl                  = no
; pbc                     = xyz  
;
gen_vel                 = yes       
gen_temp                = 300       
gen_seed                = -1
EOF

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_production.tpr
gmx mdrun -deffnm md_production