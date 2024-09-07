# ftCBS-Data
A repository for data and calculation scripts related to the manuscript:
"Removing basis set incompleteness error in finite-temperature electronic
structure calculations: Two-electron systems".

For the data points that go into the figures use the code in `plot.py` to find
the original data sources.

## Running
For running the code it is suggested to use anaconda.
The required libraries are listed in the `conda_env.yml` file.
Following the instructions:
1. `conda env create -f conda_env.yml`
2. `conda activate ftcbs`

On conda version 24.1.2 we are able to successully run `ueg.py`, `atoms.py`,
and `plot.py`.

## Files
Due to the GitHub file size limitations we are unable to provide access to all
the Hamiltonian, determinant, eigen energies, wave function, and integral
files. We have provided an example for how to generate the relevant files for
the uniform electron gas from HANDE-QMC, however the code required for this is
currently on a development branch and will be released on the public branch
during the next release cycle for HANDE-QMC. For the helium system, we have
provided the eigen energy files needed to recreate thermodynamic quantities. We
also provide the small integral files for the helium system needed by
HANDE-QMC. For access to the larger helium integral files and UEG files please
direct inquires to the corresponding author on the manuscript.

## Directories
Below we outline basic descriptions for each directory contained within this
repository.
1. atoms_data
    - Contains the analysed data used to generate the plots for the helium atom
      and hydrogen atom within the manuscript.
2. he_pbc_3x3x3
    - Contains the small integral files and HANDE-QMC output files with eigen
      energies for the smaller volume helium atom.
3. he_pbc_8x8x8
    - Contains the small integral files and HANDE-QMC output files with eigen
      energies for the larger volume helium atom.
4. rs001_M00038
    - Contains the example HANDE-QMC input files and outputs for running the
      uniform electron gas calculations to generate the files needed to
      calculate thermodynamic quantities.
5. rs001_data
    - Contains the analysed data used to generate the plots for the rs=1
      uniform electron gas system within the manuscript.
6. rs010_data
    - Contains the analysed data used to generate the plots for the rs=10
      uniform electron gas system within the manuscript.

## Citing
For referencing this data please use the manuscript which can be found at
- https://doi.org/10.1021/acs.jpca.4c03769
