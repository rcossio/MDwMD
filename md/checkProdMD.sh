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

# Restraint elimination
for file in prod.xvg mindist.xvg rmsd_first.xvg rmsd_xray.xvg gyrate.xvg sasa.xvg dssp.xvg proj.xvg
do
    scp rodper@aurora:~/$input/$file .
done

#Some controls in production MD
awk '{print $1, $2}' prod.xvg > prod_potential.xvg
../../../validation/unipd/bin/python3 ../../make_plot.py prod_potential.xvg moving_average
rm prod_potential.xvg
awk '{print $1, $3}' prod.xvg > prod_kinetic.xvg
../../../validation/unipd/bin/python3 ../../make_plot.py prod_kinetic.xvg moving_average
rm prod_kinetic.xvg
awk '{print $1, $4}' prod.xvg > prod_pressure.xvg
../../../validation/unipd/bin/python3 ../../make_plot.py prod_pressure.xvg moving_average
rm prod_pressure.xvg


#Plot PCs
../../../validation/unipd/bin/python3 ../../plot_pca.py 
cat > plot.gnu << EOF
set palette defined ( 0 "red", 1 "orange", 2 "yellow", 3 "green", 4 "blue", 5 "violet")
set size ratio -1
set grid

# Calculate the color index outside and use it directly
plot "proj.xvg" using 1:2:(\$0/1000.0) w l lc palette, "" every 100 using 1:2:(\$0*100/1000.0) w p pt 7 ps 2 lc palette notitle
pause -1
EOF
echo "Plot script created. Run with 'gnuplot plot.gnu'"