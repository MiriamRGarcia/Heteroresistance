# Heteroresistance
A **MATLAB** repository to study bacterial heteroresistance under antimicrobial stress using the model proposed in: *N. Martínez-López, C. Vilas and M. R. García (2024) A birth-death model to understand bacterial antimicrobial heteroresistance from time-kill curves.* The heteroresistant population is modelled as a multivariate Birth-Death (BD) process following the evolution of the different subpopulations.

The repository consist of the following folders and functions:

- **SSA**: This folder contains the necessary functions to generate statistically correct trajectories of the multivariate BD model proposed in the main paper using the *Stochastic Simulation Algorithm (SSA)*. The following files are included in this folder:
  - **Main_SSAmultiBD.m**: The main script where the user must define the options for simulate the multivariate BD heteroresistance model (number of trajectories, parameter values, simulation times, antimicrobial concentrations, etc.). The code permits to simulate the trajectories simultaneously for different experiments at constant antimicrobial concentrations. The results for each trajectory are saved automatically in the subfolder *Results* as a *.mat* file. Also, the user can choose to generate the trajectories using the *Direct Method (DM)* (*Gillespie, D. T. (2006) Stochastic Simulation of Chemical Kinetics*: https://doi.org/10.1146/annurev.physchem.58.032806.104637) or the *Rejection-based Stochastic Simulation Algorithm (RSSA)* (*Thanh et al. (2014) Efficient rejection-based simulation of biochemical reactions with stochastic noise and delays*: https://doi.org/10.1063/1.4896985):
  - **multiBD_SSA.m**: The function implementing DM to generate the trajectories of the multivariate BD heteroresistance model.
  - **multiBD_RSSA.m**: The function implementing the RSSA to generate the trajectories of the multivariate BD heteroresistance model.
  - **Results**: The subfolder where the simulation results are saved as *.mat* files. The user can generate the trajectories of the multivariate BD heteroresistance model from scratch or download previously generated results (used in the main paper) from the link: https://nube.iim.csic.es/s/2024-HeteroAMR-ResultsSSA. The downloaded *Results* folder must be placed in *2024-HeteroAMR/SSA*.
    
- **PE**: This folder contains the necessary functions to perform the model calibration using *Maximum Likelihood Estimation (MLE)* under the different noise scenarios described in the main paper. The following files are included in this folder:
  - **MainPE.m**: The main script where the user must define the options for MLE (type of noise, bounds on model parameters, sampling times, etc.). This function needs the SSmGO Toolbox that is available upon request (see https://www.bangalab.org/software).
  - **Functions**: The subfolder containing the necessary functions to perform MLE:
    - **costFun_MNHo.m**: Cost function for MLE in the MNHo case.
    - **costFun_MNHe.m**: Cost function for MLE in the MNHe case.
    - **costFun_PN.m**: Cost function for MLE in the PN case.
    - **odes_cte.m**: Function implementing ODEs to obtain total average counts at constant antimicrobial concentration.
    - **PlotPE.m**: Function to plot the calibration results.
      
- **OI**: This folder contains the necessary functions to perform the analysis of practical identifiability of the heteroresistance model. The following files are included in this folder:
  
- **Figures**: This folder contains the necessary functions to plot the figures included in the main paper. The following files are included in this folder:
