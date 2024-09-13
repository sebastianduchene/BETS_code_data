## Raw data and simulated data sets for: Assessing the effect of model specification and prior sensitivity on Bayesian tests of temporal signal

- cholera_aln.fasta - sequence alignment for _V. cholerae_ data set.
- powassan_alignment.fasta - sequence alignment for _Powassan virus_ data set.
- treponema.fasta - sequence alignment for _T. pallidum_ data set.
- small_simulations* - these folders contain the xml files for simulations under the corresponding prior on the population size parameter of the constant-size coalescent.
- empirical - contains a subset of xml files generated from the three sequence alignments analysed under the corresponding prior and clock model with and without sampling times.
- simulate_cholera_like contains code to generate simulations that are cholera like using R.

All analyses were conducted with BEASTv1.10.5 pre-release (now BEASTv10.5.0). No additional packages are required.
BEAGLE library 4.0.1 is also required.
As an example, the following code will assign a GPU configured by beagle-lib to run cholera_het_SC_exponential.xml .
```
/path_to_beast/bin/beast -beagle_GPU path_to_xml_file/cholera_het_SC_exponential.xml 
```
This will produce a .trees file, and .log files containing the parameter estimates.
A full guide to running BEAST can be found here: https://beast.community/getting_started
