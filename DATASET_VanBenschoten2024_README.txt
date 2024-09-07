this DATASET_VanBenschoten2024_README.txt was generated on 20241015
by WILLIAMVANBENSCHOTEN

-------------------
GENERAL INFORMATION
-------------------

Title of Dataset:
    Dataset for: Removing basis set incompleteness error in finite-temperature
    electronic structure calculations: Two-electron systems

Author Information:
    Principal Investigator: James J. Shepherd, Department of Chemistry,
    Michigan State University, East Lansing, Michigan 48824, USA,
    sheph158@msu.edu

Date of data collection:
    2024/03

Geographic location of data collection:
    University of Iowa, Chemistry Department, Iowa City, IA, USA
    Department of Chemistry, Michigan State University, East Lansing, Michigan, USA

Abstract:
    This dataset contains the data and scripts to generate the 
    figures used in the manuscript titled: "Removing basis set incompleteness
    error in finite-temperature electronic structure calculations:
    Two-electron systems". In addition to these are example calculation data
    and raw eigenenergies for the helium atom under periodic boundary
    conditions, which are used to generate the aforementioned data for the
    manuscript. For more information see the related text.

Information about funding sources or sponsorship
that supported the collection of the data:
    This research was supported by the U.S. Department of Energy,
    Office of Science, Office of Basic Energy Sciences Early Career
    Research Program (ECRP) under Award No. DE-SC0021317.
    This research also used computer resources from the University
    of Iowa and Michigan State University.

--------------------------
SHARING/ACCESS INFORMATION
--------------------------

Licenses/restrictions placed on the data, or limitations of reuse:
    GPL-3.0 license
    ODC Attribution License (ODC-By)

Recommended citation for the data:
    W. Z. Van Benschoten and J. J. Shepherd, (2024) Dataset for: Removing
    basis set incompleteness error in finite-temperature electronic structure
    calculations: Two-electron systems. (Dataset)

Citation for and links to publications that cite or use the data:
    William Z. Van Benschoten and James J. Shepherd
    The Journal of Physical Chemistry A 2024 ___ (__), ____-____
    [DOI_HERE]

Links to other publicly accessible locations of the data:
    N/A

Links/relationships to ancillary or related datasets:
    N/A

--------------------
DATA & FILE OVERVIEW
--------------------

File list (filenames, directory structure (for zipped files)
and brief description of all data files):
    Please see DATASTRUCTURE_VanBenschoten2024_README.txt
    for a full list of files and descriptions.

Relationship between files, if important for context:
    N/A

Additional related data collected that 
was not included in the current data package:
    The uniform electron gas eigen energies and functions are not
    included due to storage size limitations.

If data was derived from another source, list source:
    Uniform electron gas FCI and Hamiltonian data are generated
    with HANDE-QMC. The HANDE-QMC code can be found at:
        https://github.com/hande-qmc/hande
        (hande.uk.org) (doi:10.1021/acs.jctc.8b01217)
        Code for scanning over all symetries contained in the uniform electron
        gas will be released in a future version of HANDE-QMC.
    For the truncated analytical hydrogen atom scripts are contained within.
    For the periodic helium system integrals are generated with Vienna
    ab initio package (VASP). The VASP code can be obtained at:
        https://www.vasp.at (doi: 10.1016/0927-0256(96)00008-0)
    For more information reach out to the corresponding author on this
    dataset.

If there are there multiple versions of the dataset, list the file updated,
when and why update was made:
    N/A

--------------------------
METHODOLOGICAL INFORMATION
--------------------------

Description of methods used for collection/generation of data:
    - HANDE was used to generate and/or run the FCI and Hamiltonian files
      associated with the uniform electron gas calculations.
    - HANDE was used to run the FCI associated with the helium atom under
      periodic boundary conditions.
    - Vienna ab initio package (VASP) was used to generate integrals files
      for the helium atom under periodic boundary conditions.
    - Code for performing the sum-over-states and analytical hydrogen atoms
      is found within this repository.

Methods for processing the data:
    - Code for performing the data processing on data resulting from scripts
      within is found within the repository.

Software- or Instrument-specific information needed to interpret the data,
including software and hardware version numbers:
    - All files are readable with any modern text editor.
    - Figure generation and post-calculation analysis can be done following
      instructions in the README.md.

Standards and calibration information, if appropriate:
    N/A

Environmental/experimental conditions:
    N/A

Describe any quality-assurance procedures performed on the data:
    N/A

People involved with sample collection, processing, analysis and/or submission:
    William Z. Van Benschoten
    James J. Shepherd

--------------------------
DATA-SPECIFIC INFORMATION
--------------------------

Number of variables:
    N/A

Number of cases/rows:
    N/A

Variable list, defining any abbreviations,
units of measure, codes or symbols used:
    Please see DATADICTIONARY_VanBenschoten2024_README.txt for a full list of
    abbreviations, units of measure, codes or symbols used.
    - Unless otherwise indicated all units are in Hartree atomic units.

File Types:
    - *.pdf, figures used throughout the manuscript.
    - *.out, outputs from HANDE calculations.
    - *.csv, outputs from analysis of data.
    - *.fcidumps, files containing integrals used to define a system.
    - *.py, scripts for running python code.

Missing data codes:
    N/A

Specialized formats or other abbreviations used:
    N/A

--------------------
DATA REPRODUCIBILITY
--------------------

Calculations with HANDE-QMC:
    Code and scripts contained within this repository demonstrate how to
    reproduce the data used within the manuscript. The data needed to perform
    the sum over states calculations for the uniform electron gas are to large
    to be contained within this repository. Code for generating the needed files
    will be released at a later date within the HANDE-QMC code base.
