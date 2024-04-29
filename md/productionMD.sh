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
;nstvout                 = $2   ; velocities
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
    createMDP prod${1}.mdp
    setMDlength prod${1}.mdp $3
    setOutputFreq prod${1}.mdp 5000
    setupSimulationMode prod${1}.mdp continue
    checkFile prod${1}.mdp

    gmx grompp -f prod${1}.mdp -c $2.gro -t $2.cpt -p topol.top -o prod${1}.tpr
    checkFile prod${1}.tpr

    if [ $4 == "continue" ]; then
        gmx mdrun -deffnm prod${1} -nb gpu -cpi $2.cpt -noappend
        mv prod${1}.part000${1}.log prod${1}.log
        mv prod${1}.part000${1}.edr prod${1}.edr
        mv prod${1}.part000${1}.trr prod${1}.trr
        mv prod${1}.part000${1}.xtc prod${1}.xtc
        mv prod${1}.part000${1}.gro prod${1}.gro
    else
        gmx mdrun -deffnm prod${1} -nb gpu
    fi
    checkFiles prod${1}.log prod${1}.edr prod${1}.xtc prod${1}.gro prod${1}.cpt
}

analyseProductionMD(){

    # ENERGY & PRESSURE
    [ -f prod.xvg ] && rm prod.xvg
    for name in "$@"
    do
        gmx energy -f $name.edr -o $name.xvg -xvg none -nobackup <<< $'Pressure\nKinetic\nPotential\n0\n'
        cat $name.xvg >> prod.xvg
    done

    # SIMPLE ANALYSIS
    gmx trjcat -f "$@" -o prod.xtc -nobackup
    checkFile prod.xtc

    gmx trjconv -s prod1.tpr  -f prod.xtc -o prod_center.xtc -center -pbc mol -dt 100 -nobackup <<< $'Protein\nSystem\n'
    gmx mindist -s prod1.tpr  -f prod_center.xtc -pi -od prod_mindist.xvg -xvg none -nobackup <<< "Protein"
    gmx rms     -s nvt_f.tpr -f prod_center.xtc -o prod_rmsd_first.xvg -tu ns -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
    gmx rms     -s em.tpr    -f prod_center.xtc -o prod_rmsd_xray.xvg -tu ns -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
    gmx gyrate  -s prod1.tpr  -f prod_center.xtc -o prod_gyrate.xvg -xvg none -nobackup <<< "Protein"
    gmx sasa    -s em.tpr  -f prod_center.xtc -o prod_sasa.xvg -xvg none -nobackup <<< "Protein"

    gmx dssp    -s prod1.tpr -f prod_center.xtc  -o prod_dssp.dat -xvg none -nobackup 
    awk '{s_count = gsub(/S/, "&"); p_count = gsub(/H/, "&"); printf "%d %0.2f %0.2f\n", NR, s_count / length * 100, p_count / length * 100}' prod_dssp.dat > prod_dssp.xvg

    # PCA
    gmx covar   -s prod1.tpr -f prod_center.xtc -o prod_eigenval.xvg -v prod_eigenvect.trr -last 2 -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
    gmx anaeig  -s prod1.tpr -f prod_center.xtc -v prod_eigenvect.trr -first 1 -last 2 -2d prod_proj.xvg -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'

    for name in em nvt_r npt_r1000 npt_r200 npt_r50 npt_r10 npt_r2 npt_r1 npt_f nvt_f
    do
           gmx anaeig  -s prod1.tpr -f $name.gro -v prod_eigenvect.trr -first 1 -last 2 -2d "prod_${name}_proj.xvg" -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
    done

    # COM & MSD
    gmx trjconv -s em.tpr -f prod.xtc       -o prod_whole.xtc  -pbc whole -nobackup <<< "System"
    gmx trjconv -s em.tpr -f prod_whole.xtc -o prod_nojump.xtc -pbc nojump -nobackup <<< "System"
    rm prod_whole.xtc
    gmx traj -f prod_nojump.xtc -s em.tpr -com -ox prod_com.xvg -xvg none -nobackup <<< "Protein"
    #gmx msd -s em.tpr -f prod_nojump.xtc -o prod_msd_gmx.xvg -beginfit 0 -endfit 50000 -nobackup <<< "Protein"
    rm prod_nojump.xtc

}

checkDependency gmx
productionMD 1 nvt_f 250000000
analyseProductionMD prod1

#productionMD 2 prod1 300000000 continue
#analyseProductionMD prod1 prod2


echo "Production MD finished"
exit 0