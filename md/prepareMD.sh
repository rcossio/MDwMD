#!/bin/bash

# Useful
# MDP options: https://manual.gromacs.org/2023.2/user-guide/mdp-options.html
# File formats: https://manual.gromacs.org/2023.2/reference-manual/file-formats.html

# Details
# Force field setup: https://manual.gromacs.org/2023.2/user-guide/force-fields.html#charmm


# Advanced argument parsing
while getopts "i:h" opt; do
    case ${opt} in
        i ) input=${OPTARG} ;;
        h ) echo "Usage: prepareMD -i inputname"; exit 0 ;;
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

convertPDBToGROMACS() {
    gmx pdb2gmx -f "$input.pdb" -o protein_processed.gro -p topol.top -i posre.itp -posrefc 1000 -ff charmm36-jul2022 -water tip3p -ignh -nobackup
    cp posre.itp posre.itp.bak
    checkFiles protein_processed.gro topol.top posre.itp posre.itp.bak
}


defineSimulationBox() {
    gmx editconf -f protein_processed.gro -o protein_newbox.gro -c -d 1.5 -bt cubic -nobackup
    checkFiles protein_newbox.gro
}


solvateSystem() {
    gmx solvate -cp protein_newbox.gro -cs spc216.gro -o protein_solv.gro -p topol.top -nobackup
    checkFiles protein_solv.gro
}


addIons() {
    touch ions.mdp
    checkFile ions.mdp

    gmx grompp -f ions.mdp -c protein_solv.gro -o ions.tpr -p topol.top
    checkFile ions.tpr

    gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -conc 0.15 -neutral -nobackup <<< "SOL"
    checkFile solv_ions.gro
}


energyMinimization() {
    cat > em.mdp << EOF
; Run parameters
integrator      = steep     ; Steepest descent minimization
emtol           = 100.0     ; Tolerance for energy minimization
emstep          = 0.01      ; Minimization step size
nsteps          = 10000     ; Max steps for energy minimization

; Non-bonded parameters
cutoff-scheme   = Verlet    ; Use Verlet list for neighbor search
nstlist         = 10        ; Update neighbor list every step
rlist           = 1.2       ; Cut-off for making neighbor list (short range forces)
ns-type         = grid      ; Grid search for neighbor list
rvdw            = 1.2       ; Cut-off for Van der Waals interactions in nm

; Electrostatics
coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb        = 1.2       ; Cut-off for Coulomb interactions in nm
DispCorr        = no

; Bond constraints
constraints     = h-bonds 
constraint-algorithm = LINCS 
EOF
    checkFile em.mdp

    gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
    checkFile em.tpr

    gmx mdrun -deffnm em -nb gpu
    checkFiles em.gro em.edr em.trr em.log

    gmx energy -f em.edr -o em_potential.xvg -xvg none -nobackup <<< $'Potential\n0\n'
    checkFile em_potential.xvg
    
}

restrainedNVT() {
    createMDP nvt_r.mdp
    setMDlength nvt_r.mdp 500000
    setOutputFreq nvt_r.mdp 5000
    addPositionRestraints nvt_r.mdp
    setContinuationMode nvt_r.mdp start
    checkFile nvt_r.mdp

    gmx grompp -f nvt_r.mdp -c em.gro -r em.gro -p topol.top -o nvt_r.tpr
    checkFile nvt_r.tpr

    gmx mdrun -deffnm nvt_r -nb gpu
    checkFiles nvt_r.log nvt_r.trr nvt_r.edr nvt_r.xtc nvt_r.gro nvt_r.cpt

    gmx energy -f nvt_r.edr -o nvt_r_temperature.xvg -xvg none -nobackup <<< $'Temperature\n0\n'
    checkFile nvt_r_temperature.xvg
}

runRestrainedNPT() {
    createMDP npt_r$1.mdp
    setMDlength npt_r$1.mdp 500000
    setOutputFreq npt_r$1.mdp 5000
    addPositionRestraints npt_r$1.mdp
    addPressureCoupling npt_r$1.mdp
    setContinuationMode npt_r$1.mdp continue
    checkFile npt_r$1.mdp

    sed -e "s/1000  1000  1000/$1  $1  $1/g" posre.itp.bak > posre.itp

    gmx grompp -f npt_r$1.mdp -c $2.gro -t $2.cpt -r em.gro -p topol.top -o npt_r$1.tpr
    checkFile npt_r$1.tpr

    gmx mdrun -deffnm npt_r$1 -nb gpu
    checkFiles npt_r$1.log npt_r$1.edr npt_r$1.trr npt_r$1.xtc npt_r$1.gro npt_r$1.cpt

    gmx energy -f npt_r$1.edr -o npt_r$1.xvg -xvg none -nobackup <<< $'Pressure\nDensity\nPosition-Rest\nPotential\n0\n'
    checkFile npt_r$1.xvg
}


eliminateRestraintsInNPT() {
    runRestrainedNPT 1000 nvt_r
    runRestrainedNPT 200 npt_r1000
    runRestrainedNPT 50 npt_r200
    runRestrainedNPT 10 npt_r50
    runRestrainedNPT 2 npt_r10
}


freeNPT() {
    createMDP npt_f.mdp
    setMDlength npt_f.mdp 500000
    setOutputFreq npt_f.mdp 5000
    addPressureCoupling npt_f.mdp
    setContinuationMode npt_f.mdp continue
    checkFile npt_f.mdp

    gmx grompp -f npt_f.mdp -c npt_r2.gro -t npt_r2.cpt -p topol.top -o npt_f.tpr
    checkFile npt_f.tpr

    gmx mdrun -deffnm npt_f -nb gpu
    checkFiles npt_f.log npt_f.edr npt_f.trr npt_f.xtc npt_f.gro npt_f.cpt

    gmx energy -f npt_f.edr -o npt_f.xvg -xvg none -nobackup <<< $'Pressure\nDensity\n0\n'
    checkFile npt_f.xvg

}

freeNVT() {
    createMDP nvt_f.mdp
    setMDlength nvt_f.mdp 500000
    setOutputFreq nvt_f.mdp 5000
    setContinuationMode nvt_f.mdp continue
    checkFile nvt_f.mdp

    gmx grompp -f nvt_f.mdp -c npt_f.gro -t npt_f.cpt -p topol.top -o nvt_f.tpr
    checkFile nvt_f.tpr

    gmx mdrun -deffnm nvt_f -nb gpu
    checkFiles nvt_f.log nvt_f.edr nvt_f.trr nvt_f.xtc nvt_f.gro nvt_f.cpt

    gmx energy -f nvt_f.edr -o nvt_f.xvg -xvg none -nobackup <<< $'Pressure\n0\n'
    checkFile nvt_f.xvg

    gmx trjconv -s nvt_f.tpr  -f nvt_f.xtc -o nvt_f_center.xtc -center -pbc mol -nobackup <<< $'Protein\nSystem\n'
    checkFile nvt_f_center.xtc
    gmx mindist -s nvt_f.tpr  -f nvt_f_center.xtc -pi -od nvt_f_mindist.xvg -xvg none -nobackup <<< "Protein"
    checkFile nvt_f_mindist.xvg
}


checkDependency "gmx"

convertPDBToGROMACS
defineSimulationBox
solvateSystem
addIons
energyMinimization
restrainedNVT
eliminateRestraintsInNPT
freeNPT
freeNVT

echo "MD Preparation completed successfully."
