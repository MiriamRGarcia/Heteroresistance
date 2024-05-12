%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MainSSA_multiBD: SSA direct method to generate trajectories of the 
%                  multivariate BD heteroresistance model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables

% ----------------------------------------------------------------------- %
% User-defined settings:

% Seed to generate uniformly distributed random numbers:
seed  = 1;

% Number of SSA trajectories:
Ntraj = 1000;

% Maximum cell count to stop simulation:
N_TL  = 1e7;

% Parameter values:
b_S       = 0.63;                                                          % Natural birth rate of S cells;
b_R       = 0.36;                                                          % Natural birth rate of R cells;

alpha_b   = 2;                                                             % Shape coefficient on birth rate;

d_maxS    = 3.78;                                                          % Maximal kill rate of S cells;

alpha_d   = 3;                                                             % Shape coefficient on death rate;
beta_d    = 0.4;                                                           % Shape coefficient on death rate; 

EC_50d    = 1;                                                             % Half maximal inhibitory concentration;
H_d       = 1;                                                             % Hill coefficient;

xi_SR     = 1e-6;                                                          % Modification rate between S and R cells;
k_xi      = log(1e2);                                                      % Decay velocity on modification rate;

N_T0      = 1e6;                                                           % Initial cell count;
lambda_T0 = 50;                                                            % Decay velocity on initial heteroresistance distribution;

pars      = [b_S;b_R;alpha_b;d_maxS;alpha_d;beta_d;EC_50d;H_d;...          % Parameter array;
             xi_SR;k_xi;N_T0;lambda_T0];


% Antimicrobial concentration in each experiment (assumed constant):
MIC_S = EC_50d*(b_S/(d_maxS - b_S))^(1/H_d);                               % Minimum inhibitory concentration of S cells,
Cexp  = MIC_S*[0 1 2 4 8 32];                                              % Antimicrobial concentration,

% Time discretisation:
t0   = 0;                                                                  % Initial time [h];
tf   = 48;                                                                 % Final time [h]; 
ht   = 1e-3;                                                               % Time step;        
tmod = t0:ht:tf;                                                           % Time discretisation;    

% Discretisation of AMR level:
ra   = 0;                                                                  % Minimum AMR level (between entire sensitivity = 0 and entire resistance = 1),
rb   = 1;                                                                  % Maximum AMR level (between entire sensitivity = 0 and entire resistance = 1),
nr   = 50;                                                                 % Size of AMR level discretisation,
r    = linspace(ra, rb, nr).';                                             % AMR level discretisation,
%r    = [0.1 0.4 0.7 0.9];                                                  % Non-equispaced discretisation;


% ----------------------------------------------------------------------- %
% Call SSA and save results:
rng(seed)

for itraj = 1:Ntraj
    
    file_name = sprintf('Results/resSSA_%03u.mat', itraj);
    
    [CFUS, CFUST] = SimSSA_multiBD(tmod, r, pars, Cexp, N_TL);
    
    save(file_name, 'pars', 'tmod', 'r', 'Cexp', 'CFUS', 'CFUST', 'seed', 'N_TL')
    
end
