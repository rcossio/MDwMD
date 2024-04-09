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

checkFile() {
    [ -f $1 ] || { echo >&2 "Error: File $1 is not found.  Aborting."; exit 1; }
}

checkFiles() {
    for file in "$@"; do
        checkFile $file
    done
}

createMDP () {
cat > $1 << EOF
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

EOF
}

setMDlength () {
    cat >> $1 << EOF
; Run parameters
integrator              = md 
dt                      = 0.002     
nsteps                  = $2

EOF
}

setOutputFreq () {
    cat >> $1 << EOF
; Output control
nstxout-compressed      = $2   ; coordinates
nstvout                 = $2   ; velocities
nstenergy               = $2   ; energy saving
nstlog                  = $2   ; log file
nstcalcenergy 		    = 100    ; energie calculation

EOF
}

addPositionRestraints () {
    cat >> $1 << EOF
; Position restraints
define                  = -DPOSRES
comm-mode               = None
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

setContinuationMode () {
    if [ $2 == "start" ]; then
        cat >> $1 << EOF
; MD continuation
continuation            = no  ; Bond constraints
gen-vel                 = yes ; Velocity generation
gen-temp                = 293 ; Velocity generation
gen-seed                = -1  ; Velocity generation

EOF
    elif [ $2 == "continue" ]; then
        cat >> $1 << EOF
; MD continuation
continuation            = yes ; Bond constraints
gen-vel                 = no  ; Velocity generation

EOF
    fi
}

productionMD() {
    createMDP prod.mdp
    setMDlength prod.mdp 25000000
    setOutputFreq prod.mdp 25000
    setupSimulationMode prod.mdp continue
    checkFile prod.mdp

    gmx grompp -f prod.mdp -c nvt_f.gro -t nvt_f.cpt -p topol.top -o prod.tpr
    checkFile prod.tpr

    gmx mdrun -deffnm prod -nb gpu
    checkFiles prod.log prod.edr prod.trr prod.xtc prod.gro prod.cpt

    gmx energy -f prod.edr -o prod.xvg -xvg none -nobackup <<< $'Pressure\nKinetic\nPotential\n0\n'

    gmx trjconv -s prod.tpr  -f prod.xtc -o prod_center.xtc -center -pbc mol -nobackup <<< $'Protein\nSystem\n'
    gmx mindist -s prod.tpr  -f prod_center.xtc -pi -od mindist.xvg -xvg none -nobackup <<< "Protein"
    gmx rms     -s nvt_r.tpr -f prod_center.xtc -o rmsd_first.xvg -tu ns -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
    gmx rms     -s em.tpr    -f prod_center.xtc -o rmsd_xray.xvg -tu ns -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
    gmx gyrate  -s prod.tpr  -f prod_center.xtc -o gyrate.xvg -xvg none -nobackup <<< "Protein"
    gmx sasa    -s em.tpr  -f prod_center.xtc -o sasa.xvg -xvg none -nobackup <<< "Protein"

    gmx dssp    -s prod.tpr -f prod_center.xtc  -o dssp.dat -xvg none -nobackup 
    awk '{s_count = gsub(/S/, "&"); p_count = gsub(/H/, "&"); printf "%d %0.2f %0.2f\n", NR, s_count / length * 100, p_count / length * 100}' dssp.dat > dssp.xvg

    gmx covar   -s prod.tpr -f prod_center.xtc -o eigenval.xvg -v eigenvect.trr -last 2 -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
    gmx anaeig  -s prod.tpr -f prod_center.xtc -v eigenvect.trr -first 1 -last 2 -2d proj.xvg -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'

    for name in em nvt_r npt_r1000 npt_r200 npt_r50 npt_r10 npt_r2 npt_f nvt_f
    do
        gmx anaeig  -s prod.tpr -f $name.gro -v eigenvect.trr -first 1 -last 2 -2d "${name}_proj.xvg" -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
    done

    gmx trjconv -s em.tpr -f prod.xtc       -o prod_whole.xtc  -pbc whole -nobackup <<< "System"
    gmx trjconv -s em.tpr -f prod_whole.xtc -o prod_nojump.xtc -pbc nojump -nobackup <<< "System"
    gmx traj -f prod_nojump.xtc -s em.tpr -com -ox com.xvg -xvg none -nobackup <<< "Protein"
    gmx msd -s em.tpr -f prod_whole.xtc -o msd_gmx.xvg -beginfit 0 -endfit 50000 -nobackup <<< "Protein"

}

checkDependency gmx
productionMD


echo "Production MD finished"
exit 0



# 10. Analysis

#For calculation of diffusion:
# check: https://mattermodeling.stackexchange.com/questions/9195/how-to-calculate-diffusion-coefficient-from-msd-graph-using-gromacs
# https://gromacs.bioexcel.eu/t/diffusion-constant-using-gmx-msd/2198

# https://github.com/bio-phys/DiffusionGLS

# https://gromacs.org-gmx-users.maillist.sys.kth.narkive.com/OxMyExmH/gmx-users-unwrap-trajectory-file-using-pbc-nojump

# IMPORTANT about the COM: https://gromacs.bioexcel.eu/t/does-gmx-use-the-center-of-mass-com-to-calculate-the-msd/1968/3


#Think about this warning
# WARNING: Masses and atomic (Van der Waals) radii will be guessed

# Resources to prepare: 10min 156M(bpti) (100 frames)
# Resources to run: 1h20 1.6G(50ns bpti) trr:~1.17G xtc:~350M (5000 frames) I want to reduce it to 1000 frames
# Assuming 20 proteins: 30h, 37G (too much space, with the change is about 7.7G)