import sys
from Bio.PDB import PDBParser
from Bio.SeqUtils import IUPACData, seq1

def calculate_molecular_mass_with_adjustments(pdb_file):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    
    # Manually define a dictionary for atomic masses of common elements, including Iron (Fe)
    atomic_masses = {
        'H': 1.008,   # Hydrogen
        'D': 2.014,   # Deuterium
        'C': 12.01,   # Carbon
        'N': 14.01,   # Nitrogen
        'O': 16.00,   # Oxygen
        'P': 30.97,   # Phosphorus
        'S': 32.07,   # Sulfur
        'Fe': 55.845, # Iron
        'Se': 78.971, # Selenium
        'Zn': 65.38,  # Zinc
        'Mg': 24.305, # Magnesium
        'Ca': 40.078, # Calcium
        'Mn': 54.938, # Manganese
        'Cu': 63.546, # Copper
        'Co': 58.933, # Cobalt
	    'F': 18.998,  # Fluorine
        'Na': 22.990, # Sodium
	    'Cl': 35.45,  # Chlorine
        'K': 39.098,  # Potassium
	    'I': 126.90,  # Iodine
	    'Hg': 200.59, # Mercury
        'Cd': 112.41  # Cadmium
    }
    
    # Dictionary for additional hydrogens per residue type (simplified estimation)
    additional_hydrogens = {
        'A': 5, 'R': 11, 'N': 6, 'D': 5, 'C': 5,
        'E': 7, 'Q': 8, 'G': 3, 'H': 7, 'I': 13,
        'L': 13, 'K': 12, 'M': 9, 'F': 9, 'P': 7,
        'S': 5, 'T': 7, 'W': 10, 'Y': 9, 'V': 11
    }
    
    # Check if there are hydrogens in the PDB
    contains_hydrogens = any(atom.element == 'H' for atom in structure.get_atoms())
    
    # Initialize the total mass
    total_mass = 0
    
    # Iterate over all residues in the structure and sum their atomic masses
    for residue in structure.get_residues():
        residue_name = residue.get_resname().strip()
        try:
            one_letter_residue = seq1(residue_name)
        except KeyError:
            sys.stderr.write(f"Warning: {residue_name} is not a standard amino acid. Skipping.\n")
            continue

        for atom in residue.get_atoms():
            atom_name = atom.element.strip().capitalize()
            if atom_name in atomic_masses:
                total_mass += atomic_masses[atom_name]
            else:
                sys.stderr.write(f"Warning: Mass of {atom_name} not found. Skipping.\n")
        
        # Estimate and add the mass for missing hydrogens if no hydrogens are in the file
        if not contains_hydrogens:
            if one_letter_residue in additional_hydrogens:
                total_hydrogens_mass = additional_hydrogens[one_letter_residue] * atomic_masses['H']
                total_mass += total_hydrogens_mass

    return total_mass

pdb_file_path = sys.argv[1]
mass = calculate_molecular_mass_with_adjustments(pdb_file_path)
print(f"{mass:.1f}")

