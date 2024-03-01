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
  curl -o "${results_dir}/${pdbCode}/${pdbCode}.pdb1.gz" "https://files.rcsb.org/download/${pdbCode}.pdb1.gz"
  gzip -d "${results_dir}/${pdbCode}/${pdbCode}.pdb1.gz"
  mv "${results_dir}/${pdbCode}/${pdbCode}.pdb1" "${results_dir}/${pdbCode}/${pdbCode}.pdb"


  # Create the input.dat file with the specified template
  cat > "${results_dir}/${pdbCode}/hydropro.dat" <<EOF
${pdbCode}          !nameOfMolecule
result              !nameForOutputFile
${pdbCode}.pdb      !pdbFile
1                   !calculationType
2.9,                !AER (radius of primary elements)
-1,                 !NSIG
20.0,               !temperature (centigrade)
0.01,               !ETA (viscosity of the solvent in poises)
50000.,             !molecularWeight
0.700,              !partialSpecificVolume (cm3/g)
1.0,                !solventDensity (g/cm3)
0,                  !numberOfValuesOfQ (try0)
0,                  !numberOfIntervals (try0)
0,                  !numberOfMCTrialsForCovolume
0,                  !IDIF=1 (yes) for full diffusion tensors
*                   !End of file


EOF

  # Navigate to the created directory
  cd "${results_dir}/${pdbCode}"

  # Execute the HydroPro command
  $hydropro_command

  # Return to the original directory
  cd - > /dev/null
  diffusionCoefficient=$(grep "Translational diffusion coefficient" "${results_dir}/${pdbCode}/result-res.txt" | awk '{print $4/ 1e-8}')
  molecularWeight=$(python3 estimate_mw.py "${results_dir}/${pdbCode}/${pdbCode}.pdb")
  echo "${pdbCode},${diffusionCoefficient},${molecularWeight=$}" >> results.csv

done < "$input_file"


