%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MainPE: Main file to perform Maximum Likelihood Estimation (MLE)
%         for calibrating the average BD heteroresistance model 
%         with synthetic time-kill data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

% Add folder with necessary functions:
addpath('Functions')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Noise assumption (noise = 'MNHo'; 'MNHe'; 'PN'):
noise     = 'PN';                                                     

% If noise = 'PN', choose method used to generate the BD trajectories 
% for loading data (direct method = SSA or rejection-based = RSSA):
method    = 'SSA'; % = 'SSA'; = 'RSSA';

% Number of replicates for model calibration:
m_traj    = 3;  

% Define number of run (to not overwrite previous calibration results):
m_run     = 2;

% Name of the file with previous calibration results and synthetic data:
load_name = 'None';%'resPE_50subpop_MNHe_3traj_run1.mat';                  % Synthetic data is generated from scratch if load_name = 'None' and
                                                                           % noise = 'MNHo','MNHe'. If noise = 'PN', synthetic data must
                                                                           % have been previously generated using the MainSSA.m.

% Discretisation of the AMR level to perform calibration (equispaced):
ra   = 0;                                                                  % Minimum AMR level;
rb   = 1;                                                                  % Maximum AMR level;
m_r  = 50;                                                                 % Number of subpopulations to calibrate the model;

% Time discretisation:
t0   = 0;                                                                  % Initial simulation time;
tf   = 48;                                                                 % Final simulation time;
ht   = 1e-3;                                                               % Time step for model simulation (solve ODEs);
texp = [2 4 6 8 10 12 16 20 24 36 48];                                     % Sampling times for model calibration;
  
% Set ODEs solver (ode15s) precision:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

% ----------------------------------------------------------------------- %
% Here start the options to generate synthetic time-kill data.
% Ignore the lines within % --- % below if  
% previous calibration results are loaded from file:

% Assumptions on measurement noise (noise = 'MNHo' or 'MNHe'):
seed  = 1;                                                                 % Set seed for data reproducibility;

sd    = 0.5;                                                               % Standard deviation of the measurement error
                                                                           % in the MNHo case (log scale);
var_a = 1;                                                                 % Parameters of variance in MNHe case;
var_b = 2;

% Choose replicates to calibrate (if noise = 'PN'):
itraj = [10 50 100];

% Values of the model parameters:
b_S       = 0.63;                                                          % Birth rate of S in absence of antimicrobial;                         
b_R       = 0.36;                                                          % Birth rate of R;

alpha_b   = 2;                                                             % Shape coefficient for AMR fitness cost

d_maxS    = 3.78;                                                          % Maximal kill rate of S cells by antimicrobial;
alpha_d   = 3;                                                             % Shape coefficient of bactericide inhibition;  
beta_d    = 0.4;                                                           % Shape coefficient of bactericide inhibition;

EC_50d    = 1;                                                             % Half maximal inhibitory concentration;
H_d       = 1;                                                             % Hill coefficient;

xi_SR     = 1e-6;                                                          % Modification rate of AMR level between S and R cells;
k_xi      = log(1e2);                                                      % Shape parsameter of modification rate;

N_T0      = 1e6;                                                           % Initial population size;
lambda_T0 = 50;                                                            % Shape coefficient of initial condition;

% Discretisation of the AMR level to generate synthetic time-kill data:
m_rdata   = 50;                                                            % Number of subpopulation;

% Minimum and maximum AMR levels used to generate synthetic data
% (note that m_rdata, ra_data and rb_data has no effect if noise = 'PN',
% as data must be previously generated using MainSSA.m):
ra_data   = 0;
rb_data   = 1;

% Define array of (constant) antimicrobial concentrations:
MIC_S = EC_50d*(b_S/(d_maxS - b_S))^(1/H_d);                               % Minimum inhibitory concentration of S cells,
Cexp  = MIC_S*[0 1 2 4 8];                                                 % Array of antimicrobial concentrations,

% End of the options to generate synthetic time-kill data
% ----------------------------------------------------------------------- %

% ESS solver options:
opts.maxtime      = 60*1;                                                  % Maximum optimisation time;                                      
opts.maxeval      = 1.0e6;                                                 % Maximum number of evaluations of the cost function;
opts.strategy     = 1;                                                     % (=1) fast, (=3) robust;
opts.local.solver = 'fminsearch';                                          % Local solver;
                                                                           % |'fmincon'|'solnp'|'wdn2fb'|'fsqp';
opts.local.finish = [];                                                    % Local solver for final refinement;
opts.local.n1     = 1;                                                     % Maximum number of iterations of the local solver;
opts.local.n2     = 1;

% Set bounds on parameters for the optimisation problem:
b_Smin       = 0.5;                                                     
b_Smax       = 5.0;

b_Rmin       = 1e-2;                                                         
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

var_bmin     = 0;
var_bmax     = 2;

% Plot calibration results (plot_res = 1) or not (plot_res = 0);
plot_res = 1;                                                              

% Define colors for each antimicrobial concentration (only for plot_res=1):
col1 = [39, 183, 222]/256;                                                 
col2 = [250, 128, 114]/256;
col3 = [75, 92, 56]/256;
col4 = [250, 141, 34]/256;
col5 = [90, 14, 45]/256;
col6 = [242, 220, 35]/256;

col  = [col1;col2;col3;col4;col5;col6];

% Define markers for each antimicrobial concentration (only for plot_res=1):
mks{1} = 's';                                                              
mks{2} = 'diamond';
mks{3} = 'v';
mks{4} = '^';
mks{5} = '>';
mks{6} = '<';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user-defined settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% Name of the file to keep the calibration results:
res_name = sprintf('Results/ResPE/resPE_%usubpop_%s_%utraj_run%u.mat',...    
                    m_r, noise, m_traj, m_run);                         

%%                
% ----------------------------------------------------------------------- %
% Load calibration results or generate syntetic data for 
% MNHo and MNHe cases if not previously generated 
%(for the PN case, data must be previously generated):

load_name = sprintf('Results/ResPE/%s', load_name);

if exist(load_name, 'file')    
    
    % Load previously calibration results:
    if strcmp(noise, 'PN')
        load(load_name, 'r', 'r_data', 'tsim', 'texp', 'Cexp', 'pars', 'pars_opt', 'f_best',...
         'itraj', 'N_Tdata', 'N_Tave_data','seed', 'Weights', 'Var_data', 'ODEoptions')
    elseif strcmp(noise, 'MNHo')
        load(load_name, 'r', 'r_data', 'tsim', 'texp', 'Cexp', 'pars', 'pars_opt', 'f_best',...
            'N_Tdata','N_Tave_data','seed', 'Weights', 'Var_data', 'sd','ODEoptions')
    else
        load(load_name, 'r', 'r_data', 'tsim', 'texp', 'Cexp', 'pars', 'pars_opt', 'f_best',...
         'N_Tdata','N_Tave_data','seed', 'Weights', 'Var_data', 'var_a', 'var_b','ODEoptions')
    end
    
    fprintf('\n>> Data from file: %s was loaded sucessfully.', load_name)
    
    % Problem sizes:
    m_t     = numel(tsim);   
    m_r     = numel(r);
    m_rdata = numel(r_data);
    m_texp  = numel(texp);
    m_e     = numel(Cexp);
    
    if m_r > 2
        par_names = {'b_S','b_R','alpha_b','d_Smax','alpha_d','beta_d',...
                     'EC_50d','H_d','xi_SR', 'k_xi','N_T0','lambda_T0'};
    else
        par_names = {'b_S','b_R','d_Smax','EC_50d','H_d','xi_SR','N_T0','lambda_T0'};
    end  
    fprintf('\n>> The noise asummption is: %s', noise)
    fprintf('\n>> The number of subpopulations to perform calibration is: %u', m_r)
    fprintf('\n>> With a equispaced discretisation between r = %.2e and r = %.2e', r(1), r(m_r))
    fprintf('\n>> The number of subpopulations used to generate data is: %u', m_rdata)
    fprintf('\n>> With a equispaced discretisation between r = %.2e and r = %.2e', r_data(1), r_data(m_rdata))
    fprintf('\n>> The initial guess for the model parameters is:')   
    for ip = 1:size(par_names, 2)
        aux_par_names = cell2mat(par_names(:,ip));
        fprintf('\n>> %s ---> %.4e', aux_par_names, pars_opt(ip))
    end   
    if strcmp(noise, 'MNHe')
        fprintf('\n>> var_b ---> %.4e\n', pars_opt(end))
    end   
    fprintf('\n>> Press any button to continue with the current settings.')
    
    pause
    
    % Find indexes of sampling times within simulation times:
    texp_ind = find(ismember(tsim, texp));
           
else
    fprintf('\n>> The user have selected first run of parameter calibration.')
    fprintf('\n>> MLE will be performed without loading previous calibration results.')
    fprintf('\n>> Or maybe the load_name set by user does not exist as a valid .mat file.')
    fprintf('\n>> Syntetic time-kill data will be generated with the following settings:')    
    fprintf('\n>> The noise asummption is: %s', noise)
    fprintf('\n>> The number of subpopulations to perform calibration is: %u', m_r)
    fprintf('\n>> With a equispaced discretisation between r = %.2e and r = %.2e', ra, rb)
    fprintf('\n>> The number of subpopulations used to generate data is: %u', m_rdata)
    fprintf('\n>> With a equispaced discretisation between r = %.2e and r = %.2e', ra_data, rb_data)  
    fprintf('\n>> Press any button to continue with the current settings.')
    
    pause
    
    % Construct array of AMR levels used for calibration:                                                                      
    r = linspace(ra, rb, m_r).';
    
    % Construct array of AMR levels used for generate data:
    r_data = linspace(ra_data, rb_data, m_rdata).';
    
    % Construct simulation times used for calibration:
    tsim     = t0:ht:tf;                                                   % Array with times to simulate model;
    tsim     = sort(unique([tsim texp])).'; 
        
    texp_ind = find(ismember(tsim, texp));                                 % Indexes of experimental times within
                                                                           % simulation times;
        
    % Problem sizes:
    m_t     = numel(tsim);   
    m_texp  = numel(texp);
    m_e     = numel(Cexp);    
    
    % Load previously generated data of SSA trajectories for the PN case:
    if strcmp(noise, 'PN')
   
        traj_name = sprintf('Results/ResSSA/res%s_%03u.mat', method, itraj(1));
        
        if exist(traj_name, 'file')
            % Load AMR level, simulation times, sampling times and
            % antimicrobial concentrations:
            tsim_aux   = tsim;
            r_aux      = r;
            Cexp_aux   = Cexp;
            
            load(traj_name, 'r', 'tsim', 'Cexp', 'pars', 'N_T')
            
            tsim = tsim.';
            
            Exps = find(ismember(Cexp, Cexp_aux));
            
            if numel(tsim) ~= numel(tsim_aux) || numel(r_data) ~= numel(r) || numel(Exps) ~= numel(Cexp_aux)
                
                fprintf('\n>> The user has set simulation times, AMR levels or antimicrobial concentrations for calibration')
                fprintf('\n>> differing with those used to generate synthetic time-kill data in MainSSA.m')
                fprintf('\n>> Please, change setup for calibration and run again.')
                
                return
                
            elseif numel(find(tsim - tsim_aux)) + numel(find(r_data - r)) > 0
                fprintf('\n>> The user has set simulation times or AMR levels for calibration')
                fprintf('\n>> differing with those used to generate synthetic time-kill data in MainSSA.m')
                fprintf('\n>> Please, change setup for calibration and run again.')
                
                return          
            end
            
            r    = r_aux;
            Cexp = Cexp_aux;
            
        else
            fprintf('\n>> Calibration using PN case have been selected.')
            fprintf('\n>> However, there is no previously generated SSA data for the trajectory itraj = %u', itraj(1))
            fprintf('\n>> Calibration can not be performed. Please, generate SSA data using MainSSA.m or change trajectories to load.')
            
            return
        end
        
        % Initialise trajectories to calibrate at the sampling times:
        N_Tdata = zeros(m_texp, m_e, m_traj);
        
        N_Tdata(1:m_texp, 1:m_e, 1) = N_T(texp_ind, Exps);
        
        aux_ii = 2;
        
        for ii = itraj(2:end)
            
            traj_name = sprintf('Results/resSSA/res%s_%03u.mat', method, ii);
        
            if exist(traj_name, 'file')
                
                % Load total cell count data:
                load(traj_name, 'N_T')
                
                N_Tdata(1:m_texp, 1:m_e, aux_ii) = N_T(texp_ind, Exps);
                
            else
                fprintf('\n>> Calibration using PN case have been selected.')
                fprintf('\n>> However, there is no previously generated SSA data for the trajectory itraj = %u', ii)
                fprintf('\n>> Calibration can not be performed. Please, generate SSA data using MainSSA.m or change trajectories to load.')
            
                return
            end

            aux_ii = aux_ii + 1;
        end
        
        N_Tave_data = sum(N_Tdata, 3)/m_traj;
        
        % Define weights to remove NaN data (only for PN case):
        Weights              = ones(m_texp, m_e);
        NaN_ind              = find(isnan(N_Tave_data));
        Weights(NaN_ind)     = 0;

        N_Tave_data(NaN_ind) = 0;

        % Sample variance of replicates (only for PN case)
        Var_data = sum((N_Tdata - repmat(N_Tave_data, 1, 1, m_traj)).^2, 3)/(m_traj - 1);

        Var_data(NaN_ind) = 1;

    else
        pars = [b_S;b_R;alpha_b;d_maxS;alpha_d;beta_d;EC_50d;H_d;...       % Array with exact model parameters;
                    xi_SR;k_xi;N_T0;lambda_T0];
        % Generate random trajectories for the MNHo or MNHe cases
        % with the user-defined options:       
        R = repmat(r_data, 1, m_rdata) - repmat(r_data.', m_rdata, 1);                       
        R = R - triu(R) + tril(R).';                                       % Auxiliary matrix for ODEs;
       
        % Call to function simulating average BD model:
        if m_r < 2 || m_rdata < 2                                          % Stop calibration if m_r < 2;
            fprintf('\n>> The user has selected m_r = %u of subpopulations to calibrate.', m_r)
            fprintf('\n>> The user has selected m_rdata = %u of subpopulations to generate data.', m_rdata)
            fprintf('\n>> However, the model is not defined for m_r < 2.')
            fprintf('\n>> Please, change the number of AMR levels and run again.') 
            return     
            
        elseif m_rdata > 2   
            
            % Calculate total average counts with the current settings:
            [~, N_Tdata] = Sim_aveBD(r_data, R, tsim, Cexp, pars, ODEoptions);
        else
            if numel(find(r_data - [0;1])) + numel(find(r - [0;1])) > 0
                fprintf('\n>> The user has selected m_r = %u of subpopulations to calibrate and m_rdata = %u to generate data.', m_r, m_rdata)
                fprintf('\n>> However, the model is code is not well defined for m_r = 2 and r distinct from [0;1].')
                fprintf('\n>> Please, change the setup for calibration and run again.')
                return
            end
            
            ind_noIdentPars                = [3 5 6 10];                   % Indexes of meaningless parameters;
            ind_IdentPars                  = 1:numel(pars);
            ind_IdentPars(ind_noIdentPars) = [];
            
            % Calculate total average counts with the current settings:   
            [~, N_Tdata] = Sim_aveBD_SR(r_data, R, tsim, Cexp, pars(ind_IdentPars), ODEoptions);
        end
        
        % Add measuremennt noise to average total counts:
        rng(seed)                                                          % Set seed for generating randoms;
        
        if strcmp(noise, 'MNHo')                                           % Add homoscedastic noise to log10(average);
            N_Tdata = log10(repmat(N_Tdata, 1, 1, m_traj)) + ...
                      sd*randn(m_t, m_e, m_traj);
        else                                                               % Ad heterocedastic noise to average;
            var     = var_a*N_Tdata.^var_b;
            N_Tdata = repmat(N_Tdata, 1, 1, m_traj) + ...
                repmat(sqrt(var), 1, 1, m_traj).*randn(m_t, m_e, m_traj);     
        end  
        
        N_Tdata     = N_Tdata(texp_ind, 1:m_e, 1:m_traj);                  % Data of trajectories at sampling times;
        
        N_Tave_data = sum(N_Tdata, 3)/m_traj;                              % Average of the replicates used for calibration at the sampling times;
        
        Weights  = ones(m_texp, m_e);
        Var_data = ones(m_texp, m_e);
    end
                                             
end


%%
% ----------------------------------------------------------------------- %
% Pass options to ESS:

% Array of lower bounds for parameters necessary for ESS:
problem.x_L = [b_Smin;b_Rmin;alpha_bmin;d_maxSmin;alpha_dmin;beta_dmin;...
               EC_50dmin;H_dmin;xi_SRmin;k_ximin;N_T0min;lambda_T0min]; 

% Array of upper bounds for parameters necessary for ESS:
problem.x_U = [b_Smax;b_Rmax;alpha_bmax;d_maxSmax;alpha_dmax;beta_dmax;...
               EC_50dmax;H_dmax;xi_SRmax;k_ximax;N_T0max;lambda_T0max];

% Slight different setups when m_r = 2 and m_r > 2:
if m_r > 2
    sim_name                       = 'Sim_aveBD';                          % Name of the function simulating the model;
    pars_nom                       = pars;                                 % Nominal parameter values to scale the problem;    
else
    sim_name                       = 'Sim_aveBD_SR';                       % Name of the function simulating the model;
    ind_noIdentPars                = [3 5 6 10];                           % Indexes of meaningless parameters;
    ind_IdentPars                  = 1:numel(pars);
    ind_IdentPars(ind_noIdentPars) = [];
    problem.x_L                    =  problem.x_L(ind_IdentPars);          % Take bounds only on actual parameters;
    problem.x_U                    =  problem.x_U(ind_IdentPars);
    
    pars_nom                       = pars(ind_IdentPars);                  % Nominal parameter values to scale the problem;    
end

% Include var_b as additional parameter to calibrate for noise = 'MNHe':
if strcmp(noise, 'MNHe')
    problem.x_L = [problem.x_L;var_bmin];
    problem.x_U = [problem.x_U;var_bmax];
    
    pars_nom    = [pars_nom;var_b];
end

% Scale bounds for obtaining a better-posed calibration problem:
problem.x_L = problem.x_L./pars_nom;
problem.x_U = problem.x_U./pars_nom;


if exist(load_name, 'file')
    problem.x_0 = pars_opt./pars_nom;
else
    
    rng('shuffle')
    
    problem.x_0 = problem.x_L + (problem.x_U - problem.x_L).*rand(size(problem.x_L));
    
    
    % Print initial guess:
    x_0_auxprint = problem.x_0.*pars_nom;
    
    fprintf('\n>> The initial guess for the model parameters is:')
    if m_r > 2
        par_names = {'b_S','b_R','alpha_b','d_Smax','alpha_d','beta_d',...
                     'EC_50d','H_d','xi_SR', 'k_xi','N_T0','lambda_T0'};
    else
        par_names = {'b_S','b_R','d_Smax','EC_50d','H_d','xi_SR','N_T0','lambda_T0'};
    end 
    
    
    for ip = 1:size(par_names, 2)
        aux_par_names = cell2mat(par_names(:,ip));
        fprintf('\n>> %s ---> %.4e', aux_par_names, x_0_auxprint(ip))
    end   
    if strcmp(noise, 'MNHe')
        fprintf('\n>> var_b ---> %.4e\n', x_0_auxprint(end))
    end 

end

%%
% ----------------------------------------------------------------------- %
% Calibrate model from synthetic time-kill data using ESS:

R = repmat(r, 1, m_r) - repmat(r.', m_r, 1);                       
R = R - triu(R) + tril(R).';                                               % Auxiliary matrix for the cost function;

problem.f = sprintf('costFun_%s',noise);                                   % Name of the cost function

% Call the optimization routine (ESS):
Results  = ess_kernel(problem, opts, r, R, tsim, texp_ind, Cexp, N_Tave_data,...         
                      pars_nom, Var_data, Weights, ODEoptions, sim_name);

% Obtain optimal parameters:
pars_opt = Results.xbest.'.*pars_nom;    

% Obtain value of cost function:
f_best   = Results.fbest;                                                  

% Remove the mat file generated by ESS:
delete ess_report.mat                                                      


% ----------------------------------------------------------------------- %
% Save results:
if strcmp(noise, 'PN')
    save(res_name, 'r', 'r_data', 'tsim', 'texp', 'Cexp', 'pars', 'pars_opt', 'f_best',...
         'itraj', 'N_Tdata', 'N_Tave_data','seed', 'Weights', 'Var_data', 'ODEoptions')
else
    % Obtain estimates of additional parameters in MNHo and MNHe cases:  
    [~, N_Tmod] = feval(sim_name, r, R, tsim, Cexp, pars_opt, ODEoptions); % Total cell counts with calibrated parameters;
    
    N_Tmod      = N_Tmod(texp_ind, 1:m_e);                                 % Total cell counts at sampling times;
    
    auxN_Tdata = reshape(N_Tave_data, [], 1);
    auxN_Tmod  = reshape(N_Tmod, [], 1);
        
    if strcmp(noise, 'MNHo')
        
        % Calculate estimate of standard deviation:
        auxN_Tmod  = log10(auxN_Tmod);                                     
        var        = sum((auxN_Tdata - auxN_Tmod).^2)/(m_texp*m_e - 1);
        sd_opt     = sqrt(var);
        save(res_name, 'r', 'r_data', 'tsim', 'texp', 'Cexp', 'pars', 'pars_opt', 'f_best',...
            'N_Tdata','N_Tave_data','seed', 'Weights', 'Var_data', 'sd', 'sd_opt', 'ODEoptions')
    else
        
        % Calculate estimate of var_a:
        var_b_opt = pars_opt(end);
        var_a_opt = sum((auxN_Tdata - auxN_Tmod).^2./(auxN_Tmod.^var_b_opt))/(m_texp*m_e);
        save(res_name, 'r', 'r_data', 'tsim', 'texp', 'Cexp', 'pars', 'pars_opt', 'f_best',...
         'N_Tdata','N_Tave_data','seed', 'Weights', 'Var_data', 'var_a',...
         'var_a_opt', 'var_b', 'var_b_opt', 'ODEoptions')
    end
end

%%
% ----------------------------------------------------------------------- %
% Plot results:

if plot_res == 1
    % Transform data to log10 scale:
    if strcmp(noise, 'MNHo')
        logNT_ave_data = N_Tave_data;
    else
        logNT_ave_data = log10(N_Tave_data);
    end

    Plot_ResPE(r, R, tsim, texp, Cexp, pars_opt, logNT_ave_data, col,...
               mks, ODEoptions, sim_name)
end

% Remove paths:
rmpath('Functions')
