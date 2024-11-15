'''
Script to calculate bond delta q for all the bonds in a molecule for each
partial charge set.
	* AmberTools AM1-BCC charges
	* OpenEye AM1-BCC charges
	* OpenEye AM1-BCC ELF10
	* OpenFF NAGL charges

Instructions:
	* generate partial charges using charge_molecules.py
	* modify the variables below within the EDIT THESE VARIABLES blocks below
	* Run with conda environment created from openfe.yaml 
	
'''

import pandas as pd
from rdkit import Chem
import statistics
import csv
from pathlib import Path

# length of SEEDS list should be longer than the num_charge_sets
# if running OpenEye AM1-BCC ELF10 charges
# must use seed numbers used with charge_molecules.py 
SEEDS = [
	42, 64, 139, 28, 3, 4, 52, 35, 54, 79, 
	133, 48, 106, 109, 47, 97, 148, 19, 29, 
	92, 105, 98, 22, 27, 87, 55, 33, 69, 5,
	137, 102, 59, 1, 13, 72, 123, 122, 140, 107, 
	125, 88, 7, 30, 45, 130, 0, 86, 136, 41, 8, 
	17, 152, 51, 104, 89, 39, 96, 114, 11, 71
]
CSV_OFFSET = 2

if __name__ == "__main__":


	####################################
	#### START EDIT THESE VARIABLES ####
	####################################
	sdf_file = f"../molecules/PLB_simulation_all.sdf"
	num_charge_sets = 50
	toolkit = 'openeye'   
	charge_method = 'am1bcc'
	input_path = '/Users/megosato/Desktop/testing'
	output_path = '/Users/megosato/Desktop/testing'
	####################################
	#### END   EDIT THESE VARIABLES ####
	####################################

	all_rdmols = [mol for mol in Chem.SDMolSupplier(sdf_file, removeHs=False)]

	for rdmol in all_rdmols:
		molname = rdmol.GetProp('_Name')

		# step 1
		# load in charges csv generated from charge_molecules.py
		# create a pandas dataframe
		charge_csv = f'{input_path}/{toolkit}_{charge_method}/{toolkit}_{charge_method}_{molname}_charges.csv'
		df = pd.read_csv(charge_csv)
			

		# step 2
		# Based on the bonds in the rdkit molecule loaded by sdf determine the
		# bond dq by taking the absolute value of difference in charge at 
		# the two atoms involved in the bond
		# store the bond dq for each conformer in a dictionary by bond index
		# and the indices of the atoms involved in the bond
		bond_dq_dict_list = []
		for i in range(num_charge_sets):
		
			for j,bond in enumerate(rdmol.GetBonds()):
				bgn_idx = bond.GetBeginAtomIdx()
				end_idx = bond.GetEndAtomIdx()
			

				if i == 0:
					bond_dq_dict_list.append(dict())
					bond_dq_dict_list[j]['idx'] = bond.GetIdx()
					bond_dq_dict_list[j]['a1'] = bgn_idx
					bond_dq_dict_list[j]['a2'] = end_idx
					bond_dq_dict_list[j]['atmnum1'] = bond.GetBeginAtom().GetAtomicNum()
					bond_dq_dict_list[j]['atmnum2'] = bond.GetBeginAtom().GetAtomicNum()

				bond_dq_dict_list[j][i] = abs(df.iloc[bgn_idx,i+2] - df.iloc[end_idx,i+2])
		

		# step 4
		# generate the bond dq statistics for each bond in the molecule
		for dq_dict in bond_dq_dict_list:
			dqs = [dq_dict[i] for i in range(num_charge_sets)]
			mean = statistics.mean(dqs)
			stdev = statistics.stdev(dqs)
			var = statistics.stdev(dqs) ** 2
			dq_dict['mean'] = mean
			dq_dict['stdev'] = stdev
			dq_dict['var'] = var
			dq_dict['range'] = abs(min(dqs) - max(dqs))

		# step 5
		# save the partial charges for each conformer as well as the statistics in
		# a csv file
		csv_dir_path = Path(f'{output_path}/{toolkit}_{charge_method}_dq/')
		csv_dir_path.mkdir(exist_ok=True, parents=True)

		csv_file_path = csv_dir_path / f'{toolkit}_{charge_method}_{molname}_dq.csv'
		# Define column headers
		fieldnames = [k for k in bond_dq_dict_list[0]]
		
		# Write data to CSV file
		with open(str(csv_file_path), mode='w', newline='') as file:
			writer = csv.DictWriter(file, fieldnames=fieldnames)
			
			# Write header row
			writer.writeheader()
			
			# Write data rows
			for row in bond_dq_dict_list:
				writer.writerow(row)