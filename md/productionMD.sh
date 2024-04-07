#!/bin/bash

# Advanced argument parsing
while getopts "i:h" opt; do
    case ${opt} in
        i ) input=${OPTARG} ;;
        h ) echo "Usage: productionMD -i inputname"; exit 0 ;;
        \? ) echo "Invalid option: $OPTARG" 1>&2; exit 1 ;;
    esac
done

checkDependency() {
    command -v $1 >/dev/null 2>&1 || { echo >&2 "I require $1 but it's not installed.  Aborting."; exit 1; }
}

createMDP () {
cat > $1 << EOF
; Run parameters
integrator              = md 
dt                      = 0.002
nsteps                  = 25000000

; Output control
nstxout-compressed      = 25000   ; coordinates
nstvout                 = 25000   ; velocities
nstcalcenergy 		    = 100    ; energie calculation
nstenergy               = 25000   ; energy saving
nstlog                  = 25000   ; log file

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

; Center of mass
comm-mode               = None
EOF
}

addPositionRestraints () {
    cat >> $1 << EOF
; Position restraints
define                  = -DPOSRES
EOF
}

addPressureCoupling () {
    cat >> $1 << EOF
; Pressure coupling
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau-p                   = 2.0
ref-p                   = 1.0
compressibility         = 4.5e-5
refcoord-scaling        = com
EOF
}

setupSimulationMode () {
    if [ $2 == "start" ]; then
        cat >> $1 << EOF
; Bond constraints
continuation            = no

; Velocity generation
gen-vel                 = yes
gen-temp                = 293
gen-seed                = -1
EOF
    elif [ $2 == "continue" ]; then
        cat >> $1 << EOF
; Bond constraints
continuation            = yes
; Velocity generation
gen-vel                 = no
EOF
    fi
}

freeMD() {
    createMDP prod.mdp
    setupSimulationMode prod.mdp continue
    [ -f prod.mdp ] || exit 1

    gmx grompp -f prod.mdp -c md.gro -t md.cpt -p topol.top -o prod.tpr  -maxwarn 1 #Warning with center of mass removal
    [ -f prod.tpr ] || exit 1

    gmx mdrun -deffnm prod -nb gpu
    [ -f prod.log ] || exit 1
    [ -f prod.edr ] || exit 1
    [ -f prod.trr ] || exit 1
    [ -f prod.xtc ] || exit 1
    [ -f prod.gro ] || exit 1
    [ -f prod.cpt ] || exit 1

    gmx energy -f prod.edr -o prod.xvg -xvg none -nobackup <<< $'Pressure\nKinetic\nPotential\n0\n'

    gmx trjconv -s prod.tpr -f prod.xtc -o prod_center.xtc -center -pbc mol -nobackup <<< $'Protein\nSystem\n'
    gmx mindist -s prod.tpr -f prod_center.xtc -pi -od mindist.xvg -xvg none -nobackup <<< "Protein"
    gmx rms     -s md.tpr   -f prod_center.xtc -o rmsd_first.xvg -tu ns -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
    gmx rms     -s em.tpr   -f prod_center.xtc -o rmsd_xray.xvg -tu ns -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
    gmx gyrate  -s prod.tpr -f prod_center.xtc -o gyrate.xvg -xvg none -nobackup <<< "Protein"
    gmx sasa    -s prod.tpr -f prod_center.xtc -o sasa.xvg -xvg none -nobackup <<< "Protein"

    gmx dssp    -s prod.tpr -f prod_center.xtc  -o dssp.dat -xvg none -nobackup 
    awk '{s_count = gsub(/S/, "&"); p_count = gsub(/H/, "&"); printf "%d %0.2f %0.2f\n", NR, s_count / length * 100, p_count / length * 100}' dssp.dat > dssp.xvg

    gmx covar   -s prod.tpr -f prod_center.xtc -o eigenval.xvg -v eigenvect.trr -last 2 -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
    gmx anaeig  -s prod.tpr -f prod_center.xtc -v eigenvect.trr -first 1 -last 2 -2d proj.xvg -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'

    

}

checkDependency gmx
freeMD
echo "Production MD finished"
exit 0



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