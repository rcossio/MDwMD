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

# Bring files from the server
for file in prod.xvg mindist.xvg rmsd_first.xvg rmsd_xray.xvg gyrate.xvg sasa.xvg dssp.xvg com.xvg msd_gmx.xvg
do
    scp rodper@aurora:~/$input/$file .
done
scp rodper@aurora:~/$input/*proj.xvg .

#Energy controls
awk '{print $1, $2}' prod.xvg > prod_potential.xvg
python3 ../../make_plot.py -i prod_potential.xvg --moving_average
rm prod_potential.xvg
awk '{print $1, $3}' prod.xvg > prod_kinetic.xvg
python3 ../../make_plot.py -i prod_kinetic.xvg --moving_average
rm prod_kinetic.xvg
awk '{print $1, $4}' prod.xvg > prod_pressure.xvg
python3 ../../make_plot.py -i prod_pressure.xvg --moving_average
rm prod_pressure.xvg

#Distance control
awk '{print $1, $2}' mindist.xvg > min_selfdist.xvg
python3 ../../make_plot.py -i min_selfdist.xvg
rm min_selfdist.xvg

#RMSD control
python3 ../../make_plot.py -i rmsd_first.xvg --moving_average
python3 ../../make_plot.py -i rmsd_xray.xvg --moving_average

#Gyrate control
awk '{print $1, $2}' gyrate.xvg > radius_gyr.xvg
python3 ../../make_plot.py -i radius_gyr.xvg --moving_average
rm radius_gyr.xvg

#SASA control
python3 ../../make_plot.py -i sasa.xvg --moving_average

#DSSP control
python3 ../../make_plot.py -i dssp.xvg --cols 3

#Plot PCs
python3 ../../plot_pca.py 
cat > plot_pca.gnu << EOF
set palette defined ( 0 "red", 1 "orange", 2 "yellow", 3 "green", 4 "blue", 5 "violet")
set size ratio -1
set grid

# Calculate the color index outside and use it directly
plot "proj.xvg" using 1:2:(\$0/1000.0) w l lc palette, "" every 100 using 1:2:(\$0*100/1000.0) w p pt 7 ps 2 lc palette notitle
pause -1
EOF
echo "Plot script created for PCA. Run with 'gnuplot plot_pca.gnu'"

#Diffusion
python3 ../../calc_msd.py com.xvg > msd.xvg
cat > plot_msd.gnu << EOF
set grid
plot "msd.xvg" using 1:2 w l t "MSD", "msd_gmx.xvg" using 1:2 w l t "GMX MSD"
pause -1
EOF
echo "Plot script created for PCA. Run with 'gnuplot plot_msd.gnu'"

python3 ../../calc_diffusion.py msd.xvg > diffusion.xvg

