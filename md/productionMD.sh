#!/bin/bash

# 9. NVT production 
cat > md.mdp << EOF
; Run parameters
integrator              = md 
dt                      = 0.002     
nsteps                  = 25000000

; Output control
nstxout-compressed      = 25000  ; coordinates 1000 frames
nstvout                 = 25000  ; velocities
nstcalcenergy 		    = 100    ; energie calculation
nstenergy               = 25000  ; energy saving
nstlog                  = 25000  ; log file

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
continuation            = yes

; Velocity generation
gen-vel                 = no       

; Center of mass
comm-mode               = None
EOF
# + md.mdp

gmx grompp -f md.mdp -c rst_100.gro -t rst_100.cpt -p topol.top -o md.tpr -maxwarn 1 #Warning with center of mass removal
# + md.tpr mdout.mdp
gmx mdrun -deffnm md -nb gpu
# + md.log md.trr md.xtc md.edr md.gro md.cpt

gmx energy -f md.edr -o kinetic.xvg -xvg none -nobackup <<< $'12\n0\n'
# + md_kinetic.xvg





# 10. Analysis

#For calculation of diffusion:
# check: https://mattermodeling.stackexchange.com/questions/9195/how-to-calculate-diffusion-coefficient-from-msd-graph-using-gromacs
# https://gromacs.bioexcel.eu/t/diffusion-constant-using-gmx-msd/2198

# https://github.com/bio-phys/DiffusionGLS

# https://gromacs.org-gmx-users.maillist.sys.kth.narkive.com/OxMyExmH/gmx-users-unwrap-trajectory-file-using-pbc-nojump

gmx trjconv -f md.xtc       -s em.gro -o md_whole.xtc  -pbc whole -nobackup
gmx trjconv -f md_whole.xtc -s em.gro -o md_nojump.xtc -pbc nojump -nobackup

gmx msd -f run01.trr -s run01.tpr -s index.ndx -o msd01.xvg -tu ns -type z


#To control
# - RMSD
# - Radius of gyration
# - SASA
# - Secondary structure content (as percentages)
# - Potential energy (to check no specific high energies are present)
# - Kinetic energy (to check it is not accumulating errors)
# - PCA (to check the protein is moving aroung the minima)

#Think about this warning
# WARNING: Masses and atomic (Van der Waals) radii will be guessed

# Resources to prepare: 10min 156M(bpti) (100 frames)
# Resources to run: 1h20 1.6G(50ns bpti) trr:~1.17G xtc:~350M (5000 frames) I want to reduce it to 1000 frames
# Assuming 20 proteins: 30h, 37G (too much space, with the change is about 7.7G)