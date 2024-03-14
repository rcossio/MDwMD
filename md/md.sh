#!/bin/bash
#The script works to prepare a neutral solvated system with a cubic box:
#I couldnt make dodecahedron box work
# I am trying to figure out how to apply possitional restrainst to backbone/alpha carbons, etc
# I am trying to understand the mdp files needed for minimization, equilibration and production runs
# I am trying to verify if the steps chosen are correct (min, NPT heat, NVT unrestraining, prod)


# TO MIND:
# Ligands need to be treated differently
# grep LIG prot.pdb > lig.pdb
# grep -v LIG prot.pdb > prot_solo.pdb
# Use CGenFF to create the .str

export PATH="/opt/gromacs/2023.2/bin:$PATH"

input="5PTI"

# 1. Convert PDB to GROMACS format, selecting force field and water model
gmx pdb2gmx -f "5PTI.pdb" -o 1.protein_prepared.gro -p topology.top -ff amber99sb -water tip3p -ignh -nobackup
#Giorgia uses Charmm18 & TIP3P
# Posre itp is created here?

# 2. Define the simulation box
gmx editconf -f 1.protein_prepared.gro -o 2.box_defined.gro -c -d 1.0 -bt cubic -nobackup


# 3. Solvate the system
gmx solvate -cp 2.box_defined.gro -cs spc216.gro -o 3.system_solvated.gro -p topology.top -nobackup


# 4. Add ions to neutralize the system 
cat > ions.mdp << EOF
; ions.mdp - used as input into grompp to generate ions.tpr
integrator      = steep     ; Algorithm (steep = steepest descent minimization)
emtol           = 1000.0    ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep          = 0.01      ; Minimization step size
nsteps          = 50000     ; Maximum number of (minimization) steps to perform
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme	= Verlet    ; Buffered neighbor searching 
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
EOF

gmx grompp -f ions.mdp -c 3.system_solvated.gro -p topology.top -o 4.ions_prep.tpr 

gmx genion -s 4.ions_prep.tpr -o 4.system_neutralized.gro -p topology.top -pname NA -nname CL -neutral -nobackup << EOF
13
EOF


# 5. Creates a group for the backbone atoms (NDX and ITP)
gmx select -s 4.system_neutralized.gro -on backbone.ndx -select 'name C or name CA or name N or name O' -nobackup

# EARLY EXIT TO NOT USE GPU
exit


# 6. Energy Minimization
cat > min_params.mdp << EOF
integrator      = steep     # Steepest descent minimization
emtol           = 1000.0    # Tolerance for energy minimization
emstep          = 0.01      # Minimization step size
nsteps          = 50000     # Max steps for energy minimization
nstlist         = 1         # Update neighbor list every step
cutoff-scheme   = Verlet    # Use Verlet list for neighbor search
ns_type         = grid      # Grid search for neighbor list
coulombtype     = PME       # Particle Mesh Ewald for long-range electrostatics
rcoulomb        = 1.0       # Cut-off for Coulomb interactions in nm
rvdw            = 1.0       # Cut-off for Van der Waals interactions in nm

; Include position restraints 
define                   = -DPOSRES
posres_BACKBONE_fc       = 1000
posres_BACKBONE          = backbone.ndx
EOF
exit 

gmx grompp -f min_params.mdp -c 4.system_neutralized.gro -n backbone.ndx -p topology.top -o min.tpr
gmx mdrun -v -deffnm min -nb gpu

# 7. NVT Equilibration (with position restraints)
cat > nvt_eq.mdp << EOF
; ... (same as before) 
define                  = -DPOSRES_CA  ; Enable position restraints for alpha carbons
; Position restraint settings for alpha carbons
posres_CA_force         = 1000       ; Force constant in kJ/mol/nm^2
EOF

gmx grompp -f nvt_eq.mdp -c min.gro -p topology.top -o nvt.tpr
gmx mdrun -v -deffnm nvt -nb gpu

# 8. NPT Equilibration (also with position restrains)
cat > npt_eq.mdp <<EOF
; Run NVT equilibration to adjust temperature
integrator          = md           
nsteps              = 100000       ; ~100 ps simulation 
dt                  = 0.002        ; 2 fs timestep
constraints         = all-bonds 
cutoff-scheme       = Verlet      
ns_type             = grid         
coulombtype         = PME
rcoulomb            = 1.0
rvdw                = 1.0
tcoupl              = v-rescale    ; Temperature coupling
tc-grps             = Protein Non-Protein 
tau_t               = 0.1          
ref_t               = 298          ; Target temperature (Kelvin)
define              = -DPOSRES_CA
posres_CA_force     = 1000
EOF

gmx grompp -f npt_eq.mdp -c nvt.gro -p topology.top -o npt.tpr
gmx mdrun -v -deffnm npt -nb gpu

# 9. Gradual Restriction Release (in multiple steps)
for force_constant in 500 200 100 50 20 10; do  # Reduce force constants
   cat > release_$force_constant.mdp <<EOF 
; Run NPT equilibration to adjust pressure
integrator      = md
nsteps          = 100000       ; ~100 ps simulation 
dt              = 0.002        
constraints     = all-bonds
cutoff-scheme   = Verlet      
ns_type         = grid         
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
tcoupl          = v-rescale
tau_t           = 0.1
ref_t           = 298
pcoupl          = Berendsen    ; Or Parrinello-Rahman 
tau_p           = 2.0          ; Pressure coupling time constant 
ref_p           = 1.0          ; Reference pressure (bar)
define          = -DPOSRES_CA 
posres_CA_force = $force_constant 
EOF

   gmx grompp -f release_$force_constant.mdp -c npt.gro -p topology.top -o release_$force_constant.tpr
   gmx mdrun -v -deffnm release_$force_constant -nb gpu

   # Update starting structure for the next step
   echo 0 | gmx trjconv -s release_$force_constant.tpr -f release_$force_constant.gro -o next_start.gro -e 25000 
done

# 10. Production Run (without restraints)
cat > production_nvt.mdp << EOF
; production_nvt.mdp

; Essential Integrator settings
integrator      = md           
nsteps          = 5000000      ; ~10 ns simulation 
dt              = 0.002        ; 2 fs timestep

; Temperature Control
tcoupl          = v-rescale    ; Thermostat (example: v-rescale)
tc-grps         = Protein Non-Protein 
tau_t           = 0.1          
ref_t           = 298          ; Target temperature (Kelvin)

; Cut-offs and Electrostatics
cutoff-scheme   = Verlet      
ns_type         = grid         
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0

; Constraints
constraints     = all-bonds ; or 'h-bonds' for more flexibility

; Output Options
nstxout         = 1000          
nstvout         = 1000          
nstenergy       = 1000
nstlog          = 1000
EOF

gmx grompp -f production_nvt.mdp -c next_start.gro -p topology.top -o production.tpr
gmx mdrun -v -deffnm production -nb gpu
