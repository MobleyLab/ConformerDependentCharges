'''
Script to charge many copies of a molecule and get statistics of the partial
charge variability at each atom across all partial charge sets
	* AmberTools AM1-BCC charges
	* OpenEye AM1-BCC charges
	* OpenEye AM1-BCC ELF10
	* OpenFF NAGL charges

Instructions:
	* modify the variables below within the EDIT THESE VARIABLES blocks below
	* if running AmberTools or OpenEye AM1-BCC charges
		* Run with conda environment created from openfe.yaml 
	* if running NAGL charges
		* Run with conda environment created from nagl.yaml 
	
'''

from utils import * 
import sys
import csv
import statistics
from pathlib import Path

# length of SEEDS list should be longer than the num_charge_sets
# if running OpenEye AM1-BCC ELF10 charges
SEEDS = [
	42, 64, 139, 28, 3, 4, 52, 35, 54, 79, 
	133, 48, 106, 109, 47, 97, 148, 19, 29, 
	92, 105, 98, 22, 27, 87, 55, 33, 69, 5,
	137, 102, 59, 1, 13, 72, 123, 122, 140, 107, 
	125, 88, 7, 30, 45, 130, 0, 86, 136, 41, 8, 
	17, 152, 51, 104, 89, 39, 96, 114, 11, 71
]


if __name__ == "__main__":

	####################################
	#### START EDIT THESE VARIABLES ####
	####################################
	sdf_file = f"../molecules/PLB_simulation_all.sdf"
	# equivalent to number of conformers you want to generate
	num_charge_sets = 50
	toolkit = 'openeye'  
	charge_method = 'am1bcc'
	output_path = '/Users/megosato/Desktop/testing'

	# only applicable to am1bccelf10 charge_method
	num_confs_elf10 = 500
	####################################
	#### END   EDIT THESE VARIABLES ####
	####################################


	# Script fails if provided toolkit and charge_method are not compatible
	if not check_provided_charge_type(toolkit, charge_method):
		print("ERROR: Invalid toolkit and charge method pairing", file=sys.stderr)
		sys.exit()
	toolkit_wrapper = get_toolkit_wrapper(toolkit)

	if len(SEEDS) < num_charge_sets and charge_method == 'am1bccelf10':
		print("ERROR: More random seeds needed, must have at least as many SEEDS as num_charge_sets", file=sys.stderr)
		sys.exit()

	all_rdmols = [mol for mol in Chem.SDMolSupplier(sdf_file, removeHs=False)]

	

	for rdmol in all_rdmols:

		offmol_orig = Molecule.from_rdkit(rdmol)
		parent_outdir = Path(f"{output_path}/{offmol_orig.name}")
		parent_outdir.mkdir(exist_ok=True,parents=True)

		# step 1
		# generate conformers based off the random seed
		# save out conformers as oeb.gz
		if charge_method != 'am1bccelf10':
			rdmol_w_confs, rdmols_by_conf = get_mols_from_random_confs(rdmol, num_charge_sets, SEEDS[0])
			offmol_w_confs = Molecule.from_rdkit(rdmol_w_confs)
			oemol_w_confs = Molecule.to_openeye(offmol_w_confs)
			ofs = oechem.oemolostream(f"{str(parent_outdir)}/{offmol_orig.name}.oeb.gz")
			ofs.SetFormat(oechem.OEFormat_OEB)
			# Write the molecule to the stream
			oechem.OEWriteMolecule(ofs, oemol_w_confs)
			# Close the stream
			ofs.close()

		if charge_method == 'am1bccelf10':
			rdmols_by_conf = []
			for i in range(num_charge_sets):
				print(i)
				rdmol_w_confs_500,_ = get_mols_from_random_confs(rdmol, num_confs_elf10, SEEDS[i])
				rdmols_by_conf.append(rdmol_w_confs_500)

				offmol_w_confs = Molecule.from_rdkit(rdmol_w_confs_500)
				oemol_w_confs = Molecule.to_openeye(offmol_w_confs)
				ofs = oechem.oemolostream(f"{str(parent_outdir)}/{offmol_orig.name}_{SEEDS[i]}.oeb.gz")
				ofs.SetFormat(oechem.OEFormat_OEB)
				# Write the molecule to the stream
				oechem.OEWriteMolecule(ofs, oemol_w_confs)
				# Close the stream
				ofs.close()
				

		# step 3
		# generate the partial charges for each conformer 
		# store the partial charges for each conformer in a dictionary by atom
		# one dictionary per conformer/partial charge set
		am1bcc_partial_charge_dict_list = list()
		for i,rdmol in enumerate(rdmols_by_conf):
			offmol = Molecule.from_rdkit(rdmol)

			offmol = gen_charges_offmol(toolkit_wrapper, charge_method, offmol)

			for j,a in enumerate(offmol.atoms):
				if i == 0:
					am1bcc_partial_charge_dict_list.append(dict())
					am1bcc_partial_charge_dict_list[j]['idx'] = a.molecule_atom_index
					am1bcc_partial_charge_dict_list[j]['atmnum'] = a.atomic_number

				am1bcc_partial_charge_dict_list[j][i] = a.partial_charge.magnitude

		# step 4
		# generate the partial charge statistics for atom in the molecule
		for a_dict in am1bcc_partial_charge_dict_list:
			charges = [a_dict[i] for i in range(num_charge_sets)]
			mean = statistics.mean(charges)
			stdev = statistics.stdev(charges)
			var = statistics.stdev(charges) ** 2
			a_dict['mean'] = mean
			a_dict['stdev'] = stdev
			a_dict['var'] = var
			a_dict['range'] = abs(min(charges) - max(charges))


		# step 5
		# save the partial charges for each conformer as well as the statistics in
		# a csv file
		csv_dir_path = Path(f'{output_path}/{toolkit}_{charge_method}_charges/')
		csv_dir_path.mkdir(exist_ok=True, parents=True)

		csv_file_path = csv_dir_path / f'{toolkit}_{charge_method}_{offmol.name}_charges.csv'
		# Define column headers
		fieldnames = [k for k in am1bcc_partial_charge_dict_list[0]]
		
		# Write data to CSV file
		with open(str(csv_file_path), mode='w', newline='') as file:
			writer = csv.DictWriter(file, fieldnames=fieldnames)
			
			# Write header row
			writer.writeheader()
			
			# Write data rows
			for row in am1bcc_partial_charge_dict_list:
				writer.writerow(row)
