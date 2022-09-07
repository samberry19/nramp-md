# nramp-md
This repository contains all code to reproduce the MD analyses used in the paper "High-resolution structures map the metal import pathway in an NRAMP transporter" by Ray et al. They should also be helpful to do similar analyses for other proteins, if you're so inclined.

Code is mostly just wrappers etc. for mdtraj functions, and will require installation of mdtraj: https://www.mdtraj.org/1.9.8.dev0/index.html. It additionally relies on numpy and scikit-learn.

The scripts contained are:

* **mdtools.py**: Contains several functions used for analyzing distances, dihedrals, and water occupancies - look in it for more details.
* **rmsd.py**: Python script to calculate root mean square deviation (RMSD) over all frames of a simulation to run via command-line
* **rmsf.py**: Python script to calculate root mean square fluctuation (RMSF) over all frames of a simulation
* **water_occupancy.py**: Python script to calculate water occupancies in a site defined by a set of residues. Also has a required config csv file (an example of which is in **water_residues_config.csv**)

Other plots were made using Jupyter notebooks, which I plan to upload here eventually.

All scripts are written by [Sam Berry](https://sam-berry.com/). If you use any code from here, please cite:

[insert citation here].
