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
for name in mindist rmsd_first rmsd_xray gyrate sasa dssp com msd_gmx
do
    scp rodper@aurora:/mnt/Italia1/MDwMD/$input/prod_$name.xvg .
done
scp rodper@aurora:/mnt/Italia1/MDwMD/$input/prod.xvg .
scp rodper@aurora:/mnt/Italia1/MDwMD/$input/*proj.xvg .

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
awk '{print $1, $2}' prod_mindist.xvg > prod_min_selfdist.xvg
python3 ../../make_plot.py -i prod_min_selfdist.xvg
rm prod_min_selfdist.xvg

#RMSD control
python3 ../../make_plot.py -i prod_rmsd_first.xvg --moving_average
python3 ../../make_plot.py -i prod_rmsd_xray.xvg --moving_average

#Gyrate control
awk '{print $1, $2}' prod_gyrate.xvg > prod_radius_gyr.xvg
python3 ../../make_plot.py -i prod_radius_gyr.xvg --moving_average
rm prod_radius_gyr.xvg

#SASA control
python3 ../../make_plot.py -i prod_sasa.xvg --moving_average

#DSSP control
python3 ../../make_plot.py -i prod_dssp.xvg --cols 3

#Plot PCs
python3 ../../plot_pca.py 
cat > plot_pca.gnu << EOF
set palette defined ( 0 "red", 1 "orange", 2 "yellow", 3 "green", 4 "blue", 5 "violet")
set size ratio -1
set grid

# Calculate the color index outside and use it directly
plot "prod_proj.xvg" using 1:2:(\$0/1000.0) w l lc palette, "" every 100 using 1:2:(\$0*100/1000.0) w p pt 7 ps 2 lc palette notitle
pause -1
EOF
echo "Plot script created for PCA. Run with 'gnuplot plot_pca.gnu'"

#Diffusion
#python3 ../../calc_msd.py prod_com.xvg > prod_msd.xvg #takes too much time
cat > plot_msd.gnu << EOF
set grid
plot "prod_msd.xvg" using 1:2 w l t "MSD", "prod_msd_gmx.xvg" using 1:2 w l t "GMX MSD"
pause -1
EOF
echo "Plot script created for PCA. Run with 'gnuplot plot_msd.gnu'"

#python3 ../../calc_diffusion.py prod_msd.xvg > prod_diffusion.xvg #takes too much time
awk '{print $2, $3, $4}' prod_com.xvg > tmp.xvg
python3 ../../calc_diffusion_hummer.py 
rm tmp.xvg

