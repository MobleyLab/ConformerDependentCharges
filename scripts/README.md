# READ ME - Project Scripts

Here, we provide the scripts used to generate partial charges, assess variability, and run absolute hydration free energy calculations.

## Instructions

### Install the two conda environments needed for this project.
1. Create openfe environment from `openfe.yaml`
	* `mamba create -f openfe.yaml`
	* This environment should be used with the following scripts
		* `charge_molecules.py`
		* `calculate_bond_dq.py`
		* `prep_fe_am1bcc.py`
		* `prep_fe_elf10.py`
		* `prep_fe_nagl_pt2.py`
2. Create nagl environment from `nagl.yaml`
	* `mamba create -f nagl.yaml`
	* This environment should be used with the following scripts
		* `prep_fe_nagl_pt1.py`

### Bulk charging molecules
These are the scripts used to generate 50 random conformers of a molecule and charge them using various partial charge generation methods.
* `charge_molecules.py`
	* given `toolkit='openeye'` and `charge_method=am1bcc`
		* should be run with environment created from `openfe.yaml`
		* charges molecules using either OpenEye AM1-BCC charges or AmberTools AM1-BCC charges
		* We generate 50 random conformers
		* We charge each conformer using the chosen AM1-BCC charge toolkit
	* given `toolkit='ambertools'` and `charge_method=am1bcc`
		* should be run with environment created from `openfe.yaml`
		* charges molecules using either OpenEye AM1-BCC charges or AmberTools AM1-BCC charges
		* We generate 50 random conformers
		* We charge each conformer using the chosen AM1-BCC charge toolkit
	* given `toolkit='openeye'` and `charge_method=am1bccelf10`
		* should be run with environment created from `openfe.yaml`
		* charges molecules using OpenEye AM1-BCC ELF10 charges
		* We generate 50 sets of 500 random conformers 
			* Each set has a different random seed
		* Each set of 500 conformers is used as an input to AM1-BCC ELF10 charge generator
		* This results in 50 AM1-BCC ELF10 partial charge sets
	* given `toolkit='openff'` and `charge_method=nagl`
		* should be run with environment created from `nagl.yaml`
		* charges molecules uisng OpenForceField's NAGL charges
		* We generate 50 random conformers
			* this is in way unnecessary because NAGL charges should be conformer independent
			* however, this allows us to prove NAGL charges are conformer indpendent

### Preparation of absolute hydration free energy calculations
* AM1-BCC charges 
	* `prep_fe_am1bcc.py`
		* should be run with environment created from `openfe.yaml`
		* prepares OpenFE absolute hydration free energy calculation start files with AM1-BCC charges 
		* uses either OpenEye or AmberTools to charge molecules
		* each replicate of the calculation will be identical except in the assigned partial charges 
			* each replicate will have partial charges generated from a random conformer
			* each replicate will begin the simulation from the same 3D coordinates / conformer
* AM1-BCC ELF 10 charges 
	* `prep_fe_elf10.py`
		* should be run with environment created from `openfe.yaml`
		* prepares OpenFE absolute hydration free energy calculation start files with OpenEye AM1-BCC ELF10 charges
		* each replicate of the calculation will be identical except in the assigned partial charges 
			* each replicate will have partial charges generated from a set of 500 random conformers
			* each replicate will begin the simulation from the same 3D coordinates / conformer
* OpenFF NAGL charges
	* `prep_fe_nagl_pt1.py`
		* should be run with environment created from `nagl.yaml`
			* it is necessary to use a separate environment to generate the charges as the OpenFF toolkit version necessary to use NAGL is incompatible with the version of OpenFE used to create the simulations start files
		* generates NAGL charges and saves the charges in a `.mol2`
		* We generate 50 random conformers
			* this is in way unnecessary because NAGL charges should be conformer independent
			* however, this allows us to prove NAGL charges are conformer indpendent
	* `prep_fe_nagl_pt2.py`
		* should be run with environment created from `openfe.yaml`
		* prepares OpenFE absolute hydration free energy calculation start files with OpenFF NAGL charges
			* requires the `.mol2` files generated from `prep_fe_nagl_pt1.py`
		* each replicate of the calculation will be exactly identical including in the assigned partial charges as NAGL is a conformer independent charge generation method


### Generation of figures
* `plot_2dmol_with_qdiff.ipynb`
	* Creates a 2D image of the molecule 
	* User selects 2 partial charge sets to compare
	* Notebook highlights each atom based on the partial charge difference for that atom between the 2 partial charge sets
* `plot_2dmol_with_bonddqdiff.ipynb`
	* Creates a 2D image of the molecule 
	* User selects 2 partial charge ∆q_bond sets to compare
	* Notebook highlights each bond based on the ∆q_bond difference for that bond between the 2 partial charge ∆q_bond sets
