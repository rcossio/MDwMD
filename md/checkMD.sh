#!/bin/bash

# Advanced argument parsing
while getopts "i:h" opt; do
    case ${opt} in
        i ) input=${OPTARG} ;;
        h ) echo "Usage: checkMD -i inputname"; exit 0 ;;
        \? ) echo "Invalid option: $OPTARG" 1>&2; exit 1 ;;
    esac
done

# Exit if usage is not respected
if [ -z $input ]; then
    echo "Usage: checkMD -i inputname"; exit 1
fi


# Create a folder
[ ! -d results ] && mkdir results
cd results

[ ! -d $input ] && mkdir $input
cd $input

# Minimisation
scp rodper@aurora:~/$input/em_potential.xvg .
../../../validation/unipd/bin/python3 ../../make_plot.py em_potential.xvg 

# Restrained NVT
scp rodper@aurora:~/$input/nvt_r_temperature.xvg .
../../../validation/unipd/bin/python3 ../../make_plot.py nvt_r_temperature.xvg moving_average

# Restraint elimination in NPT
for file in npt_r1000.xvg npt_r200.xvg npt_r50.xvg npt_r10.xvg npt_r2.xvg npt_f.xvg
do
    scp rodper@aurora:~/$input/$file .
done
cat npt_r1000.xvg npt_r200.xvg npt_r50.xvg npt_r10.xvg npt_r2.xvg | grep -v -E "#|@" | awk '{print $4, $5}' > tmp1.dat
cat npt_f.xvg | grep -v -E "#|@" | awk '{print $2, $3}' > tmp2.dat

cat tmp1.dat tmp2.dat | cat -n |  awk '{print 10*$1, $2}'> npt_r_pressure.xvg
cat tmp1.dat tmp2.dat | cat -n |  awk '{print 10*$1, $3}'> npt_r_density.xvg
rm tmp1.dat tmp2.dat

../../../validation/unipd/bin/python3 ../../make_plot.py npt_r_pressure.xvg moving_average
../../../validation/unipd/bin/python3 ../../make_plot.py npt_r_density.xvg moving_average

cat npt_r1000.xvg npt_r200.xvg npt_r50.xvg npt_r10.xvg npt_r2.xvg | grep -v -E "#|@" | awk '{print $2, $3}' | cat -n |  awk '{print 10*$1, $2, $3, -$2/$3}'> restraint_full.xvg
awk '{print $1, $2}' restraint_full.xvg > npt_r_restraint.xvg
../../../validation/unipd/bin/python3 ../../make_plot.py npt_r_restraint.xvg moving_average
rm npt_r_restraint.xvg


# Controls in free NVT
scp rodper@aurora:~/$input/nvt_f.xvg .
../../../validation/unipd/bin/python3 ../../make_plot.py nvt_f.xvg moving_average


# Trajectories
for file in em.gro nvt_r.xtc npt_r1000.xtc npt_r200.xtc npt_r50.xtc npt_r10.xtc npt_r2.xtc npt_f.xtc nvt_f.xtc
do
    scp rodper@aurora:~/$input/$file .
done


# Visualisation
cat > load_visualization.vmd << EOF
# Load the initial structure file
mol new em.gro

# Load the trajectory files
mol addfile nvt_r.xtc     type xtc first 0 last -1 step 1 waitfor all
mol addfile npt_r1000.xtc     type xtc first 0 last -1 step 1 waitfor all
mol addfile npt_r200.xtc type xtc first 0 last -1 step 1 waitfor all
mol addfile npt_r50.xtc  type xtc first 0 last -1 step 1 waitfor all
mol addfile npt_r10.xtc  type xtc first 0 last -1 step 1 waitfor all
mol addfile npt_r2.xtc  type xtc first 0 last -1 step 1 waitfor all
mol addfile npt_f.xtc  type xtc first 0 last -1 step 1 waitfor all
mol addfile nvt_f.xtc  type xtc first 0 last -1 step 1 waitfor all

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