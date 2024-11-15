'''
Part TWO of the script to prepare OpenFE Absolute Solvation Free Energies with WATER 
as the solvent (in other words Absolute Hydration Free Energies) using:
	* OpenFF NAGL charges
This script takes the NAGL charges generated in PART ONE and creates OpenFE input
files with them
It is necessary to do this in 2 steps as the OpenFE version used in the paper is not
compatible with the OpenFF toolkit versions need to run NAGL charge generation.

Instructions:
	* Run with conda environment created from openfe.yaml 
	* modify the variables below within the EDIT THESE VARIABLES blocks below
'''

from utils import * 
import json

SEED = 42

if __name__ == "__main__":

	####################################
	#### START EDIT THESE VARIABLES ####
	####################################
	sdf_file = f"../../molecules/PLB_simulations_all.sdf"
	# equivalent to number of conformers you want to generate
	num_confs = 1
	num_rpts = 5
	input_path = '/Users/megosato/Desktop/testing'
	output_path = '/Users/megosato/Desktop/testing'
	####################################
	#### END   EDIT THESE VARIABLES ####
	####################################

	all_rdmols = [mol for mol in Chem.SDMolSupplier(sdf_file, removeHs=False)]

	for rdmol in all_rdmols:
		molname = rdmol.GetProp('_Name')

		# step 1
		# Load in the .mol2 file with partial charges from NAGL and the original
		# 3D coordinates
		offmol_chg = Molecule.from_file(f'{input_path}/{molname}/charged.mol2')
		print(offmol_chg.partial_charges)

		# step 2
		# Generate the small molecule component for the openFE simulation 
		smc = openfe.SmallMoleculeComponent.from_openff(offmol_chg)

		# # step 3
		# # Create openFE network for all FE clacs
		settings = get_ahfe_settings()
		network = create_network(settings, [smc])

		# step 4
		# create the directory structure to run your calculations
		parent_outdir = Path(f"{output_path}/{molname}")
		parent_outdir.mkdir(exist_ok=True,parents=True)

		transformation = [e for e in network.edges][0]
		for rpt in range(num_rpts):
			outdir = parent_outdir / str(rpt)
			outdir.mkdir(exist_ok=True)
			transformation.dump(outdir / f"ahfe.json")

		break







	