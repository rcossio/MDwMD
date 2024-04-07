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

# NVT
scp rodper@aurora:~/$input/nvt_temperature.xvg .
../../../validation/unipd/bin/python3 ../../make_plot.py nvt_temperature.xvg moving_average

# Restraint elimination
for file in npt.xvg rst_200.xvg rst_50.xvg rst_10.xvg rst_0.xvg
do
    scp rodper@aurora:~/$input/$file .
done
cat npt.xvg rst_200.xvg rst_50.xvg rst_10.xvg | grep -v -E "#|@" | awk '{print $4, $5}' > tmp1.dat
cat rst_0.xvg | grep -v -E "#|@" | awk '{print $2, $3}' > tmp2.dat

cat tmp1.dat tmp2.dat | cat -n |  awk '{print 10*$1, $2}'> pressure.xvg
cat tmp1.dat tmp2.dat | cat -n |  awk '{print 10*$1, $3}'> density.xvg
rm tmp1.dat tmp2.dat

../../../validation/unipd/bin/python3 ../../make_plot.py pressure.xvg moving_average
../../../validation/unipd/bin/python3 ../../make_plot.py density.xvg moving_average

cat npt.xvg rst_200.xvg rst_50.xvg rst_10.xvg | grep -v -E "#|@" | awk '{print $2, $3}' | cat -n |  awk '{print 10*$1, $2, $3, -$2/$3}'> restraint_full.xvg
awk '{print $1, $2}' restraint_full.xvg > restraint.xvg
../../../validation/unipd/bin/python3 ../../make_plot.py restraint.xvg moving_average
rm restraint.xvg


# Controls in free MD
scp rodper@aurora:~/$input/md.xvg .
../../../validation/unipd/bin/python3 ../../make_plot.py md.xvg moving_average


# Trajectories
for file in em.gro nvt.xtc npt.xtc rst_200.xtc rst_50.xtc rst_10.xtc rst_0.xtc md.xtc
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
mol addfile rst_0.xtc  type xtc first 0 last -1 step 1 waitfor all

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


#gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol <<< $'Protein\nProtein\n'
#gmx mindist -s md.tpr -f md_center.xtc -pi -od mindist.xvg -xvg none -nobackup <<< "Protein"
#gmx rms -s rst_0.tpr -f md_center.xtc -o rmsd_first.xvg -tu ns -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
#gmx rms -s em.tpr -f md_center.xtc -o rmsd_xray.xvg -tu ns -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
#gmx gyrate -f md_center.xtc -s md.tpr -o gyrate.xvg -xvg none -nobackup <<< "Protein"
#gmx sasa -f md_center.xtc -s md.tpr -o sasa.xvg -xvg none -nobackup <<< "Protein"

#gmx dssp -f md_center.xtc -s md.tpr -o dssp.dat -xvg none -nobackup 
#awk '{s_count = gsub(/S/, "&"); p_count = gsub(/H/, "&"); printf "%d %0.2f %0.2f\n", NR, s_count / length * 100, p_count / length * 100}' dssp.dat > dssp.xvg

#gmx covar -f md_center.xtc -s md.tpr -o eigenval.xvg -v eigenvect.trr -last 2 -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
#gmx anaeig -f md_center.xtc -s md.tpr -v eigenvect.trr -first 1 -last 2 -2d proj.xvg -xvg none -nobackup <<< $'C-alpha\nC-alpha\n'
#cat > plot.gnu << EOF
#set palette defined ( 0 "red", 1 "orange", 2 "yellow", 3 "green", 4 "blue", 5 "violet")
#plot "proj.xvg" using 1:2:(\$0/1000.0) with lines lc palette
#pause -1
#EOF
#echo "Plot script created. Run with 'gnuplot plot.gnu'"

#gmx energy -f md.edr -o md.xvg -xvg none -nobackup <<< $'Pressure\nKinetic\nPotential\n0\n'


#Some controls in production MD
#scp rodper@aurora:~/$input/md.xvg .
#awk '{print $1, $3}' md.xvg > md_pressure.xvg
#../../../validation/unipd/bin/python3 ../../make_plot.py md_pressure.xvg moving_average
#rm md_pressure.xvg
#awk '{print $1, $2}' md.xvg > md_kinetic.xvg
#../../../validation/unipd/bin/python3 ../../make_plot.py md_kinetic.xvg moving_average
#rm md_kinetic.xvg