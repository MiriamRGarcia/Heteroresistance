# Heteroresistance

A **MATLAB** repository to study bacterial heteroresistance under antimicrobial stress using the model proposed in: 

*N. Martínez-López, C. Vilas and M. R. García (2024),*

*A birth-death model to understand bacterial antimicrobial heteroresistance from time-kill curves.* 

Available at: https://dx.doi.org/10.2139/ssrn.4825402

The heteroresistant population is modelled as a multivariate Birth-Death (BD) process following the evolution of the different subpopulations.

The repository consist of the following main files:

- **MainSSA.m**: The main script where the user must define the options (number of trajectories, parameter values, simulation times, antimicrobial concentrations, etc.) for simulating stochastic trajectories of the multivariate BD heteroresistance model. The code permits to simulate the trajectories simultaneously for different experiments at constant antimicrobial concentration. The results for each trajectory are saved automatically in the folder *Results/ResSSA* as a *.mat* file. Also, the user can choose to generate the trajectories using the *Direct Method (DM)* (*Gillespie, D. T. (2006) Stochastic Simulation of Chemical Kinetics*, available at https://doi.org/10.1146/annurev.physchem.58.032806.104637) or the *Rejection-based Stochastic Simulation Algorithm (RSSA)* (*Thanh et al. (2014) Efficient rejection-based simulation of biochemical reactions with stochastic noise and delays*, available at https://doi.org/10.1063/1.4896985).

- **MainPE.m**: The main script where the user must define the options (type of noise, bounds on model parameters, sampling times, etc.) to perform calibration of the average BD heteroresistance model using Maximum Likelihood Estimation (MLE). The results for parameter calibration are saved automatically in the folder *Results/ResPE* as a *.mat* file. This function needs the SSmGO Toolbox that is available upon request (see https://www.bangalab.org/software).
  
- **MainCI.m**: The main script where the user must define the options (type of noise, confidence level, etc.) to calculate the Fisher Information Matrix (FIM) and FIM-based confidence intervals for the parameter estimates obtained with MLE. Note that, as consequence, calibration results must be previously generated using **MainPE.m**. The results of FIM and confidence intervals are saved automatically in the folder *Results/ResFIM* as a *.mat* file.

Also, the folder **Figures** contains the necessary codes to plot all the figures shown in the main paper, including simulation of the BD heteroresistance model and MLE results.

**IMPORTANT NOTE ON CODE:** The functions contained in the repository do not consider the case with only two subpopulations distinct from the entirely sensitive (S) and resistant (R) since this would require a change in the general model structure. However, the purpose of the repository is to illustrate the techniques used in the main text, not to provide a toolbox for generic use, so the codes have not been complicated to include this functionality.
