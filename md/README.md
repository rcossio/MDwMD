# Molecular Dynamics

## Prepare MD

To prepare an MD from the PDB file, we followed these steps:
- Download the PDB file from the Protein Data Bank (PDB).
- Remove HETATM (including waters)
- Remove deuterium atoms
- Set convenient title

To prepare the MD we perform the following steps:
- Generate box, add waters and ions
- Minimize the system
- NVT with restraints of 1000 kJ/mol/nm2, 1ns 
- NPT with restraints of 1000 kJ/mol/nm2, 1ns
- NPT with restraints of 200 kJ/mol/nm2, 1ns
- NPT with restraints of 50 kJ/mol/nm2, 1ns
- NPT with restraints of 10 kJ/mol/nm2, 1ns
- NPT with restraints of 2 kJ/mol/nm2, 1ns
- NPT with restraints of 0.5 kJ/mol/nm2, 1ns
- NPT without restraints, 1ns
- NVT without restraints, 1ns

To run the preparation of the MD you can follow the steps in the protocol folder that contains an SH script and MDP files. To run it all at once you can run the script `prepareMD.sh` in the terminal, in the folder where the PDB is located. To run it, you can use the following command:

```bash
bash prepareMD.sh -i <name> >prepareMD.log 2>&1 &
```
Note that if a PDB is called `1abc.pdb`, the name is `1abc`. Also note that the redirecting the output and error to prepareMD.log is optional.

To analyse the MD, we check:
- Minimization
- Temperature
- Pressure
- Density
- Visualization of the trajectory

To analyse the MD you can run:
```bash
bash checkMD.sh -i <name>
```
Keep in mind that this will create the folder ´results´ abd inside it it will put the results of the analysis (in another folder).


## Production MD
To run the production of 500ns you can run the script `productionMD.sh` in the terminal, in the folder where the PDB is located. To run it, you can use the following command:

```bash
bash productionMD.sh -i <name> >productionMD.log 2>&1 &
```
To analyse the production MD, you can run:
```bash
bash checkProductionMD.sh -i <name>
```

# Diffusivity

We used different methods to calculate the diffusion coefficient. 

## Hummer method
This method divides the trajectory in chunks and calculates the diffusivity skipping frames of the MD (it increases the skipping frame size). For every skip-frame size (and contemplating all the chunks) it calculates the diffusivity and the quality factor. To run this methods you can use the following command:

```bash
cd results/<name>
awk '{print $2,$3,$4}' prod_com.xvg > tmp.xvg
python3 ../../calc_diffusion_hummer.py
```
The outputs are called ´D_analysis.dat´ and ´D_analysis.pdf´. You can modify the values in ´calc_diffusion_hummer.py´ to obtain convergence of the diffusivity. Keep in mind that this script uses the library ´Dfit.py´ that is in the same folder.

## Covariance estimator
This method calculates the diffusivity using the covariance matrix of the displacements of the COM. To run it you can use the following command:

```bash
cd results/<name>
python3 ../../calc_diffusion_CVE.py
```

## Other methods
Scripts for other methods can be found:
- `calc_diffusion_by_avg_msd.py` for the mean square displacement method
- `calc_diffusion_from_last_frames.py` for the last frames method
- `calc_diffusion_of_time.py` for the diffusion of time method

I suggest not to use them since they have not been tested properly.

