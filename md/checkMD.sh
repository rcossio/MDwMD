#!/bin/bash


# Advanced argument parsing
while getopts "i:h" opt; do
    case ${opt} in
        i ) input=${OPTARG} ;;
        h ) echo "Usage: checkMD -i inputname"; exit 0 ;;
        \? ) echo "Invalid option: $OPTARG" 1>&2; exit 1 ;;
    esac
done


# Create a folder
mkdir $input
cd $input


# Minimisation
scp rodper@aurora:~/$input/em_potential.xvg .


# Equilibration
for file in nvt_temperature.xvg npt_pressure.xvg npt_density.xvg
do 
    scp rodper@aurora:~/$input/$file .
    ../../validation/unipd/bin/python3 ../moving_average.py $file
done


# Restraint elimination
for file in rst_1000.xvg rst_200.xvg rst_50.xvg rst_10.xvg
do
    scp rodper@aurora:~/$input/$file .
done
cat rst_1000.xvg rst_200.xvg rst_50.xvg rst_10.xvg | grep -v -E "#|@" | awk '{print $2, $3}' | cat -n |  awk '{print 20*$1, $2, $3, -$2/$3}'> rst_all.xvg


# Trajectories
for file in em.gro nvt.xtc npt.xtc rst_200.xtc rst_50.xtc rst_10.xtc
do
    scp rodper@aurora:~/$input/$file .
done


# Visualisation
cat > load_visualization.vmd << EOF
# Load the initial structure file
mol new em.gro

# Load the trajectory files
mol addfile nvt.xtc     type xtc first 0 last -1 step 1 waitfor all
mol addfile npt.xtc     type xtc first 0 last -1 step 1 waitfor all
mol addfile rst_200.xtc type xtc first 0 last -1 step 1 waitfor all
mol addfile rst_50.xtc  type xtc first 0 last -1 step 1 waitfor all
mol addfile rst_10.xtc  type xtc first 0 last -1 step 1 waitfor all

# Remove the default Lines representation
mol delrep 0 top

# Display the protein with New Cartoon representation
mol representation NewCartoon
mol color Structure
mol selection "protein"
mol addrep top

# Display water molecules (OW) with CPK representation
mol representation CPK
mol color Name
mol selection "name OW"
mol addrep top

# Display NA and CL ions with VDW representation
mol representation VDW
mol color Name
mol selection "name NA or name CL"
mol addrep top

# Set the view to orthographic
display projection Orthographic
EOF
echo "Visualisation script created. Run with 'vmd -e load_visualization.vmd'"

mv $input results/.
