'''
Part ONE of the script to prepare OpenFE Absolute Solvation Free Energies with WATER 
as the solvent (in other words Absolute Hydration Free Energies) using:
	* OpenFF NAGL charges
This script generates the NAGL charges.
It is necessary to do this in 2 steps as the OpenFE version used in the paper is not
compatible with the OpenFF toolkit versions need to run NAGL charge generation.

Instructions:
	* Run with conda environment created from nagl.yaml 
	* modify the variables below within the EDIT THESE VARIABLES blocks below
'''

from utils import * 
import sys

SEED = 42

if __name__ == "__main__":

	####################################
	#### START EDIT THESE VARIABLES ####
	####################################
	sdf_file = f"../../molecules/PLB_simulations_all.sdf"
	# equivalent to number of conformers you want to generate
	num_confs = 1
	toolkit = 'nagl'
	charge_method = 'openff-gnn-am1bcc-0.1.0-rc.1.pt'
	output_path = '/Users/megosato/Desktop/testing'
	####################################
	#### END   EDIT THESE VARIABLES ####
	####################################


	# Script fails if provided toolkit and charge_method are not compatible
	if not check_provided_charge_type(toolkit, charge_method):
		print("ERROR: Invalid toolkit and charge method pairing", file=sys.stderr)
		sys.exit()
	toolkit_wrapper = get_toolkit_wrapper(toolkit)

	all_rdmols = [mol for mol in Chem.SDMolSupplier(sdf_file, removeHs=False)]


	for rdmol in all_rdmols:

		offmol_orig = Molecule.from_rdkit(rdmol)

		# step 1
		# generate conformers based off the random seed
		rdmol_w_confs, rdmols_by_conf = get_mols_from_random_confs(rdmol, num_confs, SEED)

		# step 2
		# Convert the rdkit molecule to an openeye molecule using openff
		# so we can write out a multiconformer molecule to a .oeb.gz
		# to save for record keeping
		parent_outdir = Path(f"{output_path}/{offmol_orig.name}")
		parent_outdir.mkdir(exist_ok=True,parents=True)

		offmol_w_confs = Molecule.from_rdkit(rdmol_w_confs)
		oemol_w_confs = Molecule.to_openeye(offmol_w_confs)
		ofs = oechem.oemolostream(f"{str(parent_outdir)}/{offmol_orig.name}.oeb.gz")
		ofs.SetFormat(oechem.OEFormat_OEB)
		# Write the molecule to the stream
		oechem.OEWriteMolecule(ofs, oemol_w_confs)
		# Close the stream
		ofs.close()

		# step 3
		# Generate the small molecule component for the openFE simulation 
		smc = openfe.SmallMoleculeComponent.from_rdkit(rdmol_w_confs)

		# step 4
		# generate charges for each small molecule component
		chg_mol = gen_charges_smc(toolkit_wrapper, charge_method, smc, offmol_orig)
		offmol_chg = chg_mol.to_openff()

		# step 5 
		# save out molecule as a mol2 file with partial charges and starting 
		# 3D coordinates
		outdir = parent_outdir 
		outdir.mkdir(exist_ok=True)
		offmol_chg.to_file(f'{outdir}/charged.mol2', file_format='mol2')

		break