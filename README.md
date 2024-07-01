# Heteroresistance

A **MATLAB** repository to study bacterial heteroresistance under antimicrobial stress using the model proposed in: *N. Martínez-López, C. Vilas and M. R. García (2024) A birth-death model to understand bacterial antimicrobial heteroresistance from time-kill curves.* The heteroresistant population is modelled as a multivariate Birth-Death (BD) process following the evolution of the different subpopulations.

The repository consist of the following main files:

- **MainSSA.m**: The main script where the user must define the options (number of trajectories, parameter values, simulation times, antimicrobial concentrations, etc.) for simulating stochastic trajectories of the multivariate BD heteroresistance model. The code permits to simulate the trajectories simultaneously for different experiments at constant antimicrobial concentration. The results for each trajectory are saved automatically in the subfolder *Results/ResSSA* as a *.mat* file. Also, the user can generate the trajectories using the *Direct Method (DM)* (*Gillespie, D. T. (2006) Stochastic Simulation of Chemical Kinetics*: https://doi.org/10.1146/annurev.physchem.58.032806.104637) or the *Rejection-based Stochastic Simulation Algorithm (RSSA)* (*Thanh et al. (2014) Efficient rejection-based simulation of biochemical reactions with stochastic noise and delays*: https://doi.org/10.1063/1.4896985).

- **MainPE.m**: The main script where the user must define the options (type of noise, bounds on model parameters, sampling times, etc.) to perform calibration of the average heteroresistance model using Maximum Likelihood Estimation (MLE).  The results for parameter calibration are saved automatically in the subfolder *Results/ResPE* as a *.mat* file. This function needs the SSmGO Toolbox that is available upon request (see https://www.bangalab.org/software).
  
- **MainCI.m**: The main script where the user must define the options (type of noise, confidence level, etc.) calculating the Fisher Information Matrix (FIM) and FIM-based confidence intervals for the parameter estimates calculated with **MainPE.m**.  The results are saved automatically in the subfolder *Results/ResFIM* as a *.mat* file.

Also, the folder **Figures** contains the necessary codes to plot all the figures included in the main paper, including simulation and calibration results.

**IMPORTANT NOTE ON CODE:** The functions contained in the repository do not consider the case with only two subpopulations distinct from the entirely sensitive (S) and resistant (R) since this requires a change in the general model structure.
