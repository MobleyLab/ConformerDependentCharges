import openfe
from openfe.protocols.openmm_afe import AbsoluteSolvationProtocol
from openfe.utils import without_oechem_backend
from openff.units import unit
from openff.toolkit import Molecule
from openeye import oeomega, oechem
import gufe
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


from openff.toolkit.utils.openeye_wrapper import OpenEyeToolkitWrapper
from openff.toolkit.utils.ambertools_wrapper import AmberToolsToolkitWrapper

try:
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
    NAGL_WRAPPER_EXISTS = True
except ImportError:
    NAGL_WRAPPER_EXISTS = False

from pathlib import Path
import statistics
import csv 
import shutil
import json

def get_mols_from_random_confs(rdmol, numconfs, seed):

	rdmol_copy = Chem.Mol(rdmol)

	AllChem.EmbedMultipleConfs(rdmol_copy, numConfs=numconfs, randomSeed=seed)

	rdmols_by_conf = []

	for conf_idx in range(rdmol_copy.GetNumConformers()):


		offmol_temp = Molecule.from_rdkit(rdmol_copy)
		offmol_temp._conformers = [offmol_temp._conformers[conf_idx]]

		rdmol_temp = Molecule.to_rdkit(offmol_temp)
		rdmols_by_conf.append(rdmol_temp)

	return rdmol_copy, rdmols_by_conf

def check_provided_charge_type(toolkit, charge_method):
	valid_toolkit_charge = [
		('ambertools', 'am1bcc'), 
		('openeye', 'am1bcc'),
		('openeye', 'am1bccelf10'),
	]
	if NAGL_WRAPPER_EXISTS:
		valid_toolkit_charge += [('nagl', 'openff-gnn-am1bcc-0.1.0-rc.1.pt')]

	return (toolkit, charge_method) in valid_toolkit_charge

def get_toolkit_wrapper(toolkit):
	toolkit_wrapper_dict = {
		'ambertools': AmberToolsToolkitWrapper(),
		'openeye': OpenEyeToolkitWrapper(),
	}

	if NAGL_WRAPPER_EXISTS:
		toolkit_wrapper_dict['nagl'] = NAGLToolkitWrapper()

	return toolkit_wrapper_dict[toolkit]


def gen_charges_smc(toolkit_wrapper, charge_method, smc, offmol_orig):
	''' Generates partial charges from the conformer(s) provided by the 
		SmallMoleculeComponent (smc)
		Creates a new SmallMoleculeComponent that uses the conformer 3D coordinates from the 
		original OpenFF Molecule but with the charges from the smc
		Assumes toolkit_wrapper and charge_method have been checked to confirm their
		compatibility using check_provided_charge_type

		toolkit_wrapper: 	ToolkitWrapper object
		charge_method: 		str that can be provided to assign_partial_charges 
		smc:				SmallMoleculeComponent object containing the conformer(s) you want to charge
		offmol_orig: 		Molecule object containing the 3D conformer coordinates to be used
	'''
	offmol = smc.to_openff()
	offmol.assign_partial_charges(
		charge_method, 
		use_conformers=offmol.conformers,
		toolkit_registry=toolkit_wrapper,
	)

	offmol_orig_tmp = Molecule(offmol_orig)
	
	offmol_orig_tmp._partial_charges = offmol._partial_charges

	return openfe.SmallMoleculeComponent.from_openff(offmol_orig_tmp)

def gen_charges_offmol(toolkit_wrapper, charge_method, offmol):
	''' Generates partial charges from the conformer(s) provided by the 
		SmallMoleculeComponent (smc)
		Creates a new SmallMoleculeComponent that uses the conformer 3D coordinates from the 
		original OpenFF Molecule but with the charges from the smc
		Assumes toolkit_wrapper and charge_method have been checked to confirm their
		compatibility using check_provided_charge_type

		toolkit_wrapper: 	ToolkitWrapper object
		charge_method: 		str that can be provided to assign_partial_charges 
		offmol_orig: 		Molecule object containing the conformer to charge
	'''
	offmol.assign_partial_charges(
		charge_method, 
		use_conformers=offmol.conformers,
		toolkit_registry=toolkit_wrapper,
	)

	return offmol

def get_ahfe_settings():
	''' returns the AbsoluteSolvationProtocol settings used for our paper
	'''
	settings = AbsoluteSolvationProtocol.default_settings()
	settings.solvent_simulation_settings.equilibration_length = 500 * unit.picosecond
	settings.solvent_simulation_settings.production_length = 10000 * unit.picosecond
	settings.vacuum_simulation_settings.equilibration_length = 500 * unit.picosecond
	settings.vacuum_simulation_settings.production_length = 1000 * unit.picosecond
	if not NAGL_WRAPPER_EXISTS:
		settings.alchemical_settings.lambda_elec_windows = 7
		settings.alchemical_settings.lambda_vdw_windows = 13
		settings.alchemsampler_settings.n_repeats = 1
		settings.alchemsampler_settings.n_replicas = 20
		settings.alchemsampler_settings.online_analysis_target_error = 0.0 * unit.boltzmann_constant * unit.kelvin

	settings.vacuum_engine_settings.compute_platform = 'CUDA'
	settings.solvent_engine_settings.compute_platform = 'CUDA'

	# Set to 0 to prevent early termination due to convergence detection
	

	return settings

def create_network(settings, charged_ligands):
	''' Creates an Absolute Solvation Protocol network
	'''

	protocol = AbsoluteSolvationProtocol(settings=settings)

	stateB = openfe.ChemicalSystem({'s': openfe.SolventComponent()})

	transformations = []

	for smc in charged_ligands:
		stateA = openfe.ChemicalSystem({'l': smc, 's': openfe.SolventComponent()})   
		t = openfe.Transformation(stateA=stateA, stateB=stateB, mapping=None, protocol=protocol, name=smc.name)
		transformations.append(t)
		
	network = openfe.AlchemicalNetwork(transformations)
	return network





