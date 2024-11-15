'''
Script to prepare OpenFE Absolute Solvation Free Energies with WATER as the solvent
(in other words Absolute Hydration Free Energies) using either:
	* AmberTools AM1-BCC charges
	* OpenEye AM1-BCC charges

Instructions:
	* Run with conda environment created from openfe.yaml 
	* modify the variables below within the EDIT THESE VARIABLES blocks below
'''

from utils import * 
import sys

SEED = 42

if __name__ == "__main__":

	####################################
	#### START EDIT THESE VARIABLES ####
	####################################
	sdf_file = f"../../molecules/PLB_simulation_all.sdf"
	# equivalent to number of conformers you want to generate
	num_confs = 5 
	toolkit = 'openeye'   # can also be 'ambertools'
	charge_method = 'am1bcc'
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
		ligands = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in rdmols_by_conf]

		# step 4
		# generate charges for each small molecule component
		charged_ligands = []
		for idx,lig in enumerate(ligands):
			chg_lig = gen_charges_smc(toolkit_wrapper, charge_method, lig, offmol_orig)
			charged_ligands.append(chg_lig)

		# step 5
		# Create openFE network for all FE clacs
		settings = get_ahfe_settings()
		network = create_network(settings, charged_ligands)

		# step 6
		# create the directory structure to run your calculations
		for idx,transformation in enumerate(network.edges):
			outdir = parent_outdir / str(idx)
			outdir.mkdir(exist_ok=True)
			transformation.dump(outdir / f"ahfe.json")







	