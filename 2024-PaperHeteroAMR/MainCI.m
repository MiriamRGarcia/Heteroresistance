%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main file to calculate confidence intervals for model parameters using 
% Fisher Information Matrix (FIM) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables

addpath('Functions')
addpath('Results')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Confidence level:
confLev  = 0.9;

% Noise assumption (= 'MNHo'; 'MNHe'; 'PN');
noise    = 'PN';   

% Choose implementation to generated trajectories (direct method = SSA or rejection based = RSSA):
method   = 'RSSA'; % = 'SSA'; = 'RSSA';

% Set ODE solver precision:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);
    
% Use new defined settings (=1) or load calibration results (=0):
load_res = 1;

if load_res < 1
    % Number of replicates used in MLE to load results:
    Ntraj   = 3; 
else
    % Values of parameters to calculate FIM:
    b_S     = 0.63;                                                        % Birth rate of S in absence of antimicrobial;                         
    b_R     = 0.36;                                                        % Birth rate of R;

    alpha_b = 2;                                                           % Shape coefficient for AMR fitness cost

    d_maxS  = 3.78;                                                        % Maximal kill rate of S cells by antimicrobial;
    alpha_d = 3;                                                           % Shape coefficient of bactericide inhibition;  
    beta_d  = 0.4;                                                         % Shape coefficient of bactericide inhibition;

    EC_50d  = 1;                                                           % Half maximal inhibitory concentration;
    H_d     = 1;                                                           % Hill coefficient;

    xi_SR   = 1e-6;                                                        % Modification rate of AMR level between S and R cells;
    k_xi    = log(1e2);                                                    % Shape parsameter of modification rate;

    N_T0      = 1e6;                                                       % Initial population size;
    lambda_T0 = 50;                                                        % Shape coefficient of initial condition;

    FIM_pars  = [b_S;b_R;alpha_b;d_maxS;alpha_d;beta_d;EC_50d;H_d;...
                 xi_SR;k_xi;N_T0;lambda_T0];
    if strcmp(noise, 'MNHo')
        sd = 0.5;
    elseif strcmp(noise, 'MNHe')
        var_a = 1;
        var_b = 2;
    else
        itraj = [1 100 1000];
    end
    
    % Time discretisation:
    t0   = 0;                                                              % Initial simulation time;
    tf   = 48;                                                             % Final simulation time;
    ht   = 1e-3;                                                           % Time step for model simulation (solve ODEs);
    tsim = t0:ht:tf;                                                       % Times to simulate ODEs;
    texp = [2 4 6 8 10 12 16 20 24 36 48];                                 % Sampling times for model calibration;


    % Discretisation of AMR level:
    equi = 1;                                                              % Equispaced discretisation of AMR level (=1) or not (=0)
                                                                           % If equi = 0 the user must define the AMR level discretisation
                                                                           % as a column array. For example: r = [0.3;0.35;0.7;0.99];
    r    = [];                                                             % Discretisation of AMR level if equi = 0;

    if equi == 1 
        ra   = 0;                                                          % Minimum AMR level;
        rb   = 1;                                                          % Maximum AMR level;
        nr   = 50;                                                         % Number of subpopulations;
        hr   = 1/(nr - 1);
        r    = (ra:hr:rb).';
    elseif numel(r) == 0
        fprintf('\n >> Please, define a discretisation of the AMR level or set "equi" option equal to one.\n')
        return
    end
    
    % Drug concentration in each experiment:
    MIC_S = EC_50d*(b_S/(d_maxS - b_S))^(1/H_d);
    Cexp  = MIC_S*[0.0 1.0 2.0 4.0 8.0].';
   
end

% Nominal parameter values for better conditioning of FIM:
pars_nom  = [0.1;0.1;1;1;1;0.1;1;1;1e-6;1;1e6;10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user-defined settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if load_res < 1
    % Load optimal parameters and calibration settings:
    results_name = sprintf('../PE/Results/resPE_%s_%utraj.mat', noise, Ntraj); 
    if strcmp(noise, 'MNHo')
        load(results_name, 'r', 'tsim', 'Cexp', 'texp', 'pars_opt', 'Weights', 'Var_data', 'sd')
    elseif strcmp(noise, 'MNHe')
        load(results_name, 'r', 'tsim', 'Cexp', 'texp', 'pars_opt', 'Weights', 'Var_data', 'var_a', 'var_b')
        pars_opt = pars_opt(1:end-1);
    else
        load(results_name, 'r', 'tsim', 'Cexp', 'texp', 'pars_opt', 'Weights', 'Var_data')
    end
    FIM_pars  = pars_opt;
    N_T0      = FIM_pars(end-1);
    lambda_T0 = FIM_pars(end);
else
    tsim     = sort(unique([tsim texp])).';   
end

texp_ind = find(ismember(tsim, texp));

% Problem sizes:
nr    = numel(r);
nt    = numel(tsim);    
ntexp = numel(texp);
np    = numel(FIM_pars);
Nexp  = numel(Cexp);

Weights = ones(ntexp, Nexp);
NaN_ind = find(Weights == 0);

% ----------------------------------------------------------------------- %
% Calculate sensitivity matrix:

% Initial condition:
f0  = exp(-lambda_T0*r);
f0  = f0/sum(f0);
N_0 = N_T0*f0;

% Call function to calculate sensitivities:
[N_T, sensN_T] = SensMultiExp(tsim, r, Cexp, FIM_pars, N_0, ODEoptions);

% ----------------------------------------------------------------------- %
% Build covariance matrix:
if strcmp(noise, 'MNHo')
    CovMatrix = sd^2*ones(ntexp*Nexp,ntexp*Nexp);
    pars_var  = [];
elseif strcmp(noise, 'MNHe')
    auxN_T    = reshape(N_T(texp_ind, 1:Nexp).', [], 1); 
    CovMatrix = diag(var_a*auxN_T.^var_b);
    pars_var  = [var_a;var_b];
else
    % Calculate sample variances:
    Ntraj   = numel(itraj);
    N_Tdata = zeros(ntexp, Nexp, Ntraj);
    aux_ii  = 1;
    
    for ii = itraj
        traj_name = sprintf('../SSA/Results/res%s_%03u.mat', method, ii);
        load(traj_name, 'N_T')
        N_Tdata(1:ntexp, 1:Nexp, aux_ii) = N_T(texp_ind, 1:Nexp);
        aux_ii = aux_ii + 1;
    end
    
    N_Tave_data = sum(N_Tdata, 3)/Ntraj;
    NaN_ind     = find(isnan(N_Tave_data));
    
    Var_data    = sum((N_Tdata - repmat(N_Tave_data, 1, 1, Ntraj)).^2, 3)/(Ntraj - 1);
    
    N_Tave_data(NaN_ind) = [];
    Var_data(NaN_ind)    = [];
    
    Var_data  = reshape(Var_data, [], 1);
    CovMatrix = diag(Var_data);
    pars_var  = [];
end

% ----------------------------------------------------------------------- %
% Build sensitivity matrix:
SensMatrix = [];

if strcmp(noise, 'MNHo')
    % Log10 of total average count:
    y_T = log10(N_T);

    % Calculate sensitivities of y_T to parameters:
    sensy_T = zeros(nt, np, Nexp);

    for ip = 1:np
        aux = reshape(sensN_T(1:nt, ip, 1:Nexp), nt, Nexp);
        aux = (aux./N_T)/log(10);
        sensy_T(1:nt, ip, 1:Nexp) = aux;
        
        % Obtain sensitivities at the sampling times:
        aux = reshape(sensy_T(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    
        % Scaling by the nominal parameter (transformation par* = par/nom_par):
        aux = pars_nom(ip)*aux;
    
        SensMatrix = [SensMatrix reshape(aux.', [], 1)];
    end
else
    for ip = 1:np
    
        % Obtain sensitivities at the sampling times:
        aux = reshape(sensN_T(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    
        aux(NaN_ind) = [];

        % Scaling by the nominal parameter (transformation par* = par/nom_par):
        aux = pars_nom(ip)*aux;
    
        SensMatrix = [SensMatrix reshape(aux.', [], 1)];
    end
end

% Normalised sensitivity matrix using the sampling times:      
% Important! this is the normalised sensitivity with respect parameters
% before reescaling using nominal values. 
normSensMatrix = [];

if strcmp(noise, 'MNHo')
     for ip = 1:np
        aux = reshape(sensy_T(texp_ind, ip, 1:Nexp), ntexp, Nexp);
        aux = FIM_pars(ip)*aux./N_T(texp_ind, 1:Nexp);
        normSensMatrix = [normSensMatrix reshape(aux.',[],1)];
    end   
else
    for ip = 1:np
        aux = reshape(sensN_T(texp_ind, ip, 1:Nexp), ntexp, Nexp);
        aux = FIM_pars(ip)*aux./N_T(texp_ind, 1:Nexp);
        normSensMatrix = [normSensMatrix reshape(aux.',[],1)];
    end
end


% ----------------------------------------------------------------------- %
% Calculate confidence intervals from FIM:
[FIM, FIMConfInt] = FIM_ConfInt(FIM_pars, pars_nom, pars_var, SensMatrix, normSensMatrix, CovMatrix, confLev, noise, Nexp, ntexp);  

% ----------------------------------------------------------------------- %
% Save results:
res_name = sprintf('Results/resFIM_%s.mat', noise);
save(res_name, 'r', 'tsim', 'texp', 'Cexp', 'FIM_pars', 'pars_nom', 'CovMatrix', 'FIM', 'FIMConfInt')

rmpath('Functions')
rmpath('Results')
