#!/bin/bash

# Define variables
results_dir="results"
hydropro_command="../../hydropro10"
input_file="pdb_list.txt"

# Ensure the results directory exists
mkdir -p "$results_dir"

echo "pdbCode,diffusionCoefficient,estimatedMolecularWeight" > results.csv
# Read each line from pdblist.txt
while IFS= read -r pdbCode; do
  # Create a directory for the pdbCode within the results directory
  mkdir -p "${results_dir}/${pdbCode}"
  cd "${results_dir}/${pdbCode}"
  
  # Download the PDB file
  curl -o $pdbCode.pdb1.gz "https://files.rcsb.org/download/${pdbCode}.pdb1.gz"
  gzip -d $pdbCode.pdb1.gz
  mv $pdbCode.pdb1 $pdbCode.pdb

  # Remove hydrogen, deuterium, water and other non-atom lines from the PDB file
  # HYDROPRO counts HETATM for the calculation
  grep -v "          H" $pdbCode.pdb | grep -v "          D" | grep -v " HOH" | grep -v " DOD" | grep -E "MODEL|ATOM|TER|HETATM|ENDMDL|END" > "${pdbCode}_clean.pdb"

  # Create the input.dat file with the specified template
  cat > hydropro.dat <<EOF
${pdbCode}          !nameOfMolecule
result              !nameForOutputFile
${pdbCode}_clean.pdb  !pdbFile
1                   !calculationType
2.9,                !AER (radius of primary elements)
-1,                 !NSIG
20.0,               !temperature (centigrade)
0.01,               !ETA (viscosity of the solvent in poises)
50000.,             !dummy molecularWeight
0.700,              !dummy partialSpecificVolume (cm3/g)
1.0,                !solventDensity (g/cm3)
0,                  !numberOfValuesOfQ (try0)
0,                  !numberOfIntervals (try0)
0,                  !numberOfMCTrialsForCovolume
0,                  !IDIF=1 (yes) for full diffusion tensors
*                   !End of file


EOF

  # Execute the HydroPro command
  $hydropro_command

  # Return to the original directory
  cd - > /dev/null
  diffusionCoefficient=$(grep "Translational diffusion coefficient" "${results_dir}/${pdbCode}/result-res.txt" | awk '{print $4/ 1e-8}')
  molecularWeight=$(python3 estimate_mw.py "${results_dir}/${pdbCode}/${pdbCode}.pdb")
  echo "${pdbCode},${diffusionCoefficient},${molecularWeight=$}" >> results.csv

done < "$input_file"


