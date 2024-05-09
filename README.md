# Heteroresistance
A MATLAB repository to study bacterial heteroresistance using the model proposed in: <i N. Martínez-López, C. Vilas and M. R. García (2024) A birth-death model of antimicrobial heteresistance in bacteria and identifiability analysis from time-kill curves.> Bacterial heteroresistance is modelled as a multivariate Birth-Death (BD) process following the evolution of the different subpopulations.

The repository consist of the following folders and functions:
- <i>SSA<i>: This folder contains the necessary functions to generate exact trajectories of the multivariate BD model proposed in the main paper using the <i>Stochastic Simulation Algorithm (SSA)<i>. Note that the provided code implement the SSA direct method, which is the less computationally efficient, as the objetive is not to develop fast simulations. The following files are included in this folder:
  - <i>MainSSA_multiBD.m<i>: The main script where the user must define the options for simulate the multivariate BD heteroresistance model (number of trajectories, parameter values, simulation times, antimicrobial concentrations, etc.). The code permits to simulate the trajectories simultaneously for different constant antimicrobial concentrations. The results are saved automatically in the folder <i>Results<i> as a <i>.mat<i> file.
  - <i>SimSSA_multiBD.m<i>: The function implementing the SSA direct method to generate the trajectories of the multivariate BD heteroresistance model.
  - Results: The folder where the simulation results are saved as <i>.mat<i> files. The user can generate the trajectories of the multivariate BD heteroresistance model from scratch or download previously generated results, used in the main paper, from the link: https://nube.iim.csic.es/s/2024-HeteroAMR-ResultsSSA. The downloaded <i>Results<i> folder must be placed in <i>2024-HeteroAMR/SSA<i>.
- <i>PE<i>: This folder contains the necessary functions to perform the model calibration using <i>Maximum Likelihood Estimation (MLE)<i> under the different noise scenarios described in the main paper. 
  - <i>MainPE.m<i>: The main script where the user must define the options for MLE (type of noise, bounds on model parameters, sampling times, etc.). This function needs the SSmGO Toolbox that is available upon request (see https://www.bangalab.org/software).
