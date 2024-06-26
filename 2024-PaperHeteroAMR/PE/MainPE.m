%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MainPE: Fit the heteroresistance model to synthetic time-kill data using 
% Maximum Likelihood Estimation (MLE) with Enhaced Scatter Search (ESS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all
clc

addpath('Functions')
addpath('Results')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Options to generate random data:
noise    = 'MNHo';                                                         % Noise assumption (= 'MNHo'; 'MNHe'; 'PN');

% Define number of run (to not overwrite previous results):
m_run    = 1;
 
% Choose method used to generated trajectories for PN case
% (direct method = SSA or rejection based = RSSA):
method   = 'RSSA'; % = 'SSA'; = 'RSSA';

m_traj   = 3;                                                              % Number of replicates;
load_res = 0;                                                              % Load previous fit results (=1) or not (=0);

if strcmp(noise, 'PN')
    itraj = [10 50 100];                                                   % Trajectories of the BD process to load;
elseif strcmp(noise, 'MNHo')
    sd = 0.5;                                                              % Standard deviation of the measurement error (log scale);
else
    var_a = 1;                                                             % Parameters of variance for MNHe;
    var_b = 2;
end

% Values of model parameters (to generate data in MNHo and MNHe cases):
b_S     = 0.63;                                                            % Birth rate of S in absence of antimicrobial;                         
b_R     = 0.36;                                                            % Birth rate of R;

alpha_b = 2;                                                               % Shape coefficient for AMR fitness cost

d_maxS  = 3.78;                                                            % Maximal kill rate of S cells by antimicrobial;
alpha_d = 3;                                                               % Shape coefficient of bactericide inhibition;  
beta_d  = 0.4;                                                             % Shape coefficient of bactericide inhibition;

EC_50d  = 1;                                                               % Half maximal inhibitory concentration;
H_d     = 1;                                                               % Hill coefficient;

xi_SR   = 1e-6;                                                            % Modification rate of AMR level between S and R cells;
k_xi    = log(1e2);                                                        % Shape parsameter of modification rate;

N_T0      = 1e6;                                                           % Initial population size;
lambda_T0 = 50;                                                            % Shape coefficient of initial condition;

seed      = 10;                                                            % Seed used to generate random data for MNHo and MNHe cases;

% Time discretisation:
t0   = 0;                                                                  % Initial simulation time;
tf   = 48;                                                                 % Final simulation time;
ht   = 1e-3;                                                               % Time step for model simulation (solve ODEs);
texp = [2 4 6 8 10 12 16 20 24 36 48];                                     % Sampling times for model calibration;


% Discretisation of AMR level (equispaced):
ra   = 0;                                                                  % Minimum AMR level;
rb   = 1;                                                                  % Maximum AMR level;
m_r  = 2;                                                                  % Number of subpopulations;

% Define array of (constant) antimicrobial concentrations:
MIC_S = EC_50d*(b_S/(d_maxS - b_S))^(1/H_d);                               % Minimum inhibitory concentration of S cells,
Cexp  = MIC_S*[0 1 2 4 8];                                                 % Array of antimicrobial concentrations,

% Set ODE solver precision:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

% ESS solver options:
opts.maxtime      = 1.0e2;                                                 % Maximum optimisation time;                                      
opts.maxeval      = 1.0e6;                                                 % Maximum number of evaluations of the cost function;
opts.strategy     = 1;                                                     % (=1) fast, (=3) robust;
opts.local.solver = 'fmincon'; %'fmincon';'solnp';  %'wdn2fb';  %'fsqp';   % Local solver;
opts.local.finish = [];%'fminsearch';%'fmincon'; %'wdn2fb'; %'fsqp';       % Local solver for final refinment;
%opts.local.n1 = 1;                                                        % Maximum number of iterations of the local solver;
%opts.local.n2 = 1;

% Plot results:
plot_res = 1;                                                              % Plot calibration results (=1) or not (=0);

% Define colors for each antimicrobial concentration:
col1 = [39, 183, 222]/256;                                                 
col2 = [250, 128, 114]/256;
col3 = [75, 92, 56]/256;
col4 = [250, 141, 34]/256;
col5 = [90, 14, 45]/256;
col6 = [242, 220, 35]/256;

col  = [col1;col2;col3;col4;col5;col6];

% Define markers for each antimicrobial concentration:
mks{1} = 's';                                                              
mks{2} = 'diamond';
mks{3} = 'v';
mks{4} = '^';
mks{5} = '>';
mks{6} = '<';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user-defined settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results_name = sprintf('Results/resPE_%usubpop_%s_%utraj_run%u.mat',...    % Name of the file to keep the PE results;
                        m_r, noise, m_traj, m_run);                         

pars = [b_S;b_R;alpha_b;d_maxS;alpha_d;beta_d;EC_50d;H_d;xi_SR;...         % Array with model parameters;
       k_xi;N_T0;lambda_T0];

tsim     = t0:ht:tf;                                                       % Array with times to simulate model;
tsim     = sort(unique([tsim texp])).'; 
texp_ind = find(ismember(tsim, texp));  

r  = linspace(ra, rb, m_r).';                                              % Array of AMR levels;


R = repmat(r, 1, m_r) - repmat(r.', m_r, 1);                               % Auxiliary matrix to solve ODEs;
R = R - triu(R) + tril(R).';

m_t    = numel(tsim);                                                      % Problem sizes;
m_texp = numel(texp);
m_e    = numel(Cexp);

% ----------------------------------------------------------------------- %
% Generate syntetic data for MNHo and MNHe if not previously generated
% (for PN case, data must be previously generated using SSA)

if load_res == 1
    
    load(results_name, 'N_Tdata', 'N_Tave_data', 'pars_opt', 'Weights', 'Var_data', 'seed');
    
else
    
    % ------------------------------------------------------------------- %
    % Load data of SSA trajectories for the PN case:
    if strcmp(noise, 'PN')
        N_Tdata = zeros(m_texp, m_e, m_traj);
        aux_ii = 1;
        for ii = itraj
            traj_name = sprintf('../SSA/Results/res%s_%03u.mat', method, ii);
            load(traj_name, 'N_T')
            N_Tdata(1:m_texp, 1:m_e, aux_ii) = N_T(texp_ind, 1:m_e);
            aux_ii = aux_ii + 1;
        end
        
    else
        
    % ------------------------------------------------------------------- %
    % Generate random trajectories for the MNHo or MNHe cases:
    
        rng(seed)                                                          % Set seed for generate randoms (MNHo and MNHe cases);
        
        N_T = zeros(m_t, m_e);                                             % Initialice total population size;        
        
        f_0 = exp(-lambda_T0*r);
        f_0 = f_0/sum(f_0);
        N_0 = N_T0*f_0;                                                    % Initial condition for ODEs;
                                                 
        Xi = xi_SR*exp(k_xi*(1 - R));                                      % Matrix with modification rates in AMR level;
        Xi = Xi - diag(diag(Xi));

        AA_aux = Xi.' - diag(sum(Xi, 2));                                  % Initialice coefficient matrix;  
        
        b     = b_S*b_R./(b_R + r.^alpha_b*(b_S - b_R));                   % Birth rate;
        d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./...
                (beta_d^alpha_d + r.^alpha_d);                             % Maximal kill rate;
        
        for iexp = 1:m_e
    
            C  = Cexp(iexp);  
            HC = C^H_d/(C^H_d + EC_50d^H_d);
            d  = d_max*HC;                                                 % Death rate;

            AA = AA_aux + diag(b - d);                                     % Coefficient matrix for ODEs;

            [~, xout] = ode15s(@(t,s) Odes_cte(t, s, AA), tsim, N_0, ODEoptions);
   
            
            N_T(1:m_t, iexp) = sum(xout, 2);                                 % Total population size;
        end

        % Add noise to data:
        N_T_aux = N_T(texp_ind, 1:m_e);
        if strcmp(noise, 'MNHo')
            N_Tdata = log10(repmat(N_T_aux, 1, 1, m_traj)) + sd*randn(m_texp, m_e, m_traj);
        else
            var     = var_a*N_T_aux.^var_b;
            N_Tdata = repmat(N_T_aux, 1, 1, m_traj) + repmat(sqrt(var), 1, 1, m_traj).*randn(m_texp, m_e, m_traj);     
        end   
    end
    
    N_Tave_data = sum(N_Tdata, 3)/m_traj;
    
    % Define weights to remove NaN data (for PN case):
    Weights = ones(m_texp, m_e);
    NaN_ind = find(isnan(N_Tave_data));
    Weights(NaN_ind) = 0;
    
    % Sample variance of replicates (for PN case)
    Var_data = sum((N_Tdata - repmat(N_Tave_data, 1, 1, m_traj)).^2, 3)/(m_traj - 1);
end

% ----------------------------------------------------------------------- %
% Calibrate model with ESS:
b_Smin       = 0.5;                                                        % Bounds for the (common) decision variables;
b_Smax       = 5.0;

b_Rmin       = 1e-1;                                                         
b_Rmax       = b_Smin - 1e-2;

alpha_bmin   = 0.1;
alpha_bmax   = 10;

d_maxSmin    = 1;                                                           
d_maxSmax    = 10;

beta_dmin    = 0.1;
beta_dmax    = 1;

alpha_dmin   = 0.1;                                                        
alpha_dmax   = 10;

EC_50dmin    = 1*MIC_S;                                                      
EC_50dmax    = 10*MIC_S;
    
H_dmin       = 0.1;                                                         
H_dmax       = 10;

xi_SRmin     = 1e-7;                                                           
xi_SRmax     = 1e-4;

k_ximin      = 1;
k_ximax      = 20;

N_T0min      = 1e4;
N_T0max      = 1e7;

lambda_T0min = 10;
lambda_T0max = 100;

problem.x_L = [b_Smin;b_Rmin;alpha_bmin;d_maxSmin;alpha_dmin;beta_dmin;EC_50dmin;H_dmin;xi_SRmin;k_ximin;N_T0min;lambda_T0min]; 
problem.x_U = [b_Smax;b_Rmax;alpha_bmax;d_maxSmax;alpha_dmax;beta_dmax;EC_50dmax;H_dmax;xi_SRmax;k_ximax;N_T0max;lambda_T0max];

if m_r < 3
    ind_noIdentPars                = [3 5 6 10];                           % Indexes of non-identifiable parameters;
    ind_IdentPars                  = 1:numel(pars);
    ind_IdentPars(ind_noIdentPars) = [];
    problem.x_L =  problem.x_L(ind_IdentPars);                             % Take bounds only on identifiable parameters;
    problem.x_U =  problem.x_U(ind_IdentPars);
end

if strcmp(noise, 'MNHe')
    var_bmin    = 0;
    var_bmax    = 2;
    problem.x_L = [problem.x_L;var_bmin];
    problem.x_U = [problem.x_U;var_bmax];
end

pars_nom    = problem.x_L + (problem.x_U - problem.x_L).*rand(size(problem.x_L));
problem.x_L = problem.x_L./pars_nom;
problem.x_U = problem.x_U./pars_nom;
np          = numel(problem.x_L);

if load_res == 1                                                           % Initial guess;
    problem.x_0 = pars_opt./pars_nom;
else
    rng('shuffle')
    problem.x_0 = problem.x_L + (problem.x_U - problem.x_L).*rand(size(problem.x_L));
end

if m_r > 2
    problem.f = sprintf('costFun_%s',noise);                               % Name of the cost function; 
else
    problem.f = sprintf('costFun_%usubpop_%s', m_r, noise);
end

Results = ess_kernel(problem, opts, r, R, tsim, texp_ind, Cexp,...         % Call the optimization function (ESS):
                     N_Tave_data, pars_nom, Var_data, Weights, ODEoptions);

pars_opt = Results.xbest.'.*pars_nom;                                      % Obtain optimal parameters;
f_best   = Results.fbest;                                                  % Obtain value of cost function;
                
delete ess_report.mat                                                      % Remove the mat file generated by ESS

% ----------------------------------------------------------------------- %
% Save results:
if strcmp(noise, 'PN')
    save(results_name, 'r', 'tsim', 'texp', 'Cexp', 'pars', 'pars_opt', 'f_best',...
         'itraj', 'N_Tdata', 'N_Tave_data','seed', 'Weights', 'Var_data', 'ODEoptions')
else
    % Calculate total counts with optimal parameters:
    if m_r > 2
        b_S       = pars_opt(1);                                                      
        b_R       = pars_opt(2);                                                     

        alpha_b   = pars_opt(3);                                                         

        d_maxS    = pars_opt(4);                                                    

        alpha_d   = pars_opt(5);                                                          
        beta_d    = pars_opt(6);                                                        

        EC_50d    = pars_opt(7);                                                        
        H_d       = pars_opt(8);                                                           

        xi_SR     = pars_opt(9);                                                      
        k_xi      = pars_opt(10);                                                   
        N_T0      = pars_opt(11);                                                       
        lambda_T0 = pars_opt(12);   
    
   
        f_0 = exp(-lambda_T0*r);
        f_0 = f_0/sum(f_0);
        N_0 = N_T0*f_0;                                                 
   
        N_Tmod = zeros(m_texp, m_e);                                     
                                                            
        Xi = xi_SR*exp(k_xi*(1 - R));
        Xi = Xi - diag(diag(Xi));

        AA_aux = Xi' - diag(sum(Xi, 2));                                  

        b     = b_S*b_R./(b_R + r.^alpha_b*(b_S - b_R));                  
        d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./...
            (beta_d^alpha_d + r.^alpha_d);   
    else
        b_S       = pars_opt(1);                                                      
        b_R       = pars_opt(2);                                                     

        d_maxS    = pars_opt(3);                                                    

        EC_50d    = pars_opt(4);                                                        
        H_d       = pars_opt(5);                                                           

        xi_SR     = pars_opt(6);                                                                                                      
        N_T0      = pars_opt(7);                                                       
        lambda_T0 = pars_opt(8);   
    
   
        f_0 = exp(-lambda_T0*r);
        f_0 = f_0/sum(f_0);
        N_0 = N_T0*f_0;                                                 
   
        N_Tmod = zeros(m_texp, m_e);                                     
                                                            
        Xi     = [0 xi_SR;xi_SR 0];

        AA_aux = Xi.' - diag(sum(Xi, 2));                                  

        b     = [b_S;b_R];                
        d_max = [d_maxS;0];
    end
    
    for iexp = 1:m_e
        C  = Cexp(iexp);  
        HC = C^H_d/(C^H_d + EC_50d^H_d);
        d  = d_max*HC;                                             

        AA = AA_aux + diag(b - d);                                    

        [~, xout] = ode15s(@(t,s) Odes_cte(t, s, AA), tsim, N_0, ODEoptions);
        
        N_Tmod(1:m_texp, iexp) = sum(xout(texp_ind,1:m_r), 2);            
        
    end    
    
    auxN_Tdata = reshape(N_Tave_data, [], 1);
    auxN_Tmod  = reshape(N_Tmod, [], 1);
        
    if strcmp(noise, 'MNHo')
        % Calculate estimate of standard deviation:
        auxN_Tmod  = log10(auxN_Tmod);
        var        = sum((auxN_Tdata - auxN_Tmod).^2)/(m_texp*m_e - 1);
        sd         = sqrt(var);
        save(results_name, 'r', 'tsim', 'texp', 'Cexp', 'pars', 'pars_opt', 'f_best',...
            'N_Tdata','N_Tave_data','seed', 'Weights', 'Var_data', 'sd','ODEoptions')
    else
        % Calculate estimate of parameter var_a:
        var_b = pars_opt(end);
        var_a = sum((auxN_Tdata - auxN_Tmod).^2./(auxN_Tmod.^var_b))/(m_texp*m_e);
        save(results_name, 'r', 'tsim', 'texp', 'Cexp', 'pars', 'pars_opt', 'f_best',...
         'N_Tdata','N_Tave_data','seed', 'Weights', 'Var_data', 'var_a', 'var_b','ODEoptions')
    end
end

% ----------------------------------------------------------------------- %
% Plot results:

if plot_res == 1
    % Transform to log10 scale:
    if strcmp(noise, 'MNHo')
        logNT_ave_data = N_Tave_data;
    else
        logNT_ave_data = log10(N_Tave_data);
    end

    PlotPE(tsim, texp, r, Cexp, pars_opt, logNT_ave_data, col, mks, ODEoptions)
end

rmpath('Functions')
rmpath('Results')
