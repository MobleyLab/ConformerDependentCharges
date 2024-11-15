# READ ME - Simulation Start Files

Here, we have gathered the simulation start files to run the absolute hydration free energy calculations from our paper. 

* The simulation start files were created using OpenFreeEnergy's `openfe v0.15.0` 
	* As this version of `openfe` is a pre-stable release version, some of the keyword arguments and setting names may change in the future, and thus our files may need to be modified
	* Each simulation is run using `openfe`'s `AbsoluteSolvationProtocol`, an absolute solvation free energy calculation protocol, using water as the solvent. 
	* `openfe.yaml` is an example of the conda environment packages used for our file preparation

* Each directory is named using the following convention 
	* `{MoleculeDatabase}_{ChargeType}`
		* `MoleculeDatabase`: either PLB (ProteinLigandBenchmarks) or FreeSolv
		* `ChargeType`: one of am1bcc, nagl, elf10
	* Lastly we append `_m1` or `_elf10` for our hardware based calculations


* Within each directory above are the names of all molecules simulated with that `MoleculeDatabase` and `ChargeType` (and hardware)
	* Within each molecule name directory are 5 directories number `0/`, `1/`, `2/`, `3/`, `4/`
		* There is one exception to this in the FreeSolv_am1bcc molecules.
		* Represents a replicate of the same molecule, where the only thing differing in the start files is the assigned partial charges.
		* `ahfe.json`: contains all the necessary information to run our `openfe` simulations
		* `run_openfe.sh`: an example SLURM script to run the `openfe` calculation using `openfe quickrun`