%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit the heteroresistance model to synthetic time-kill assay using 
% Enhaced Scatter Search (ESS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

addpath('Functions')
addpath('Results')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined options:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% Options to generate random data:
noise    = 'MNHo';                                                         % Noise assumption (= 'MNHo'; 'MNHe'; 'PN');
Ntraj    = 3;                                                              % Number of trajectories to generate data;
load_res = 0;                                                              % Load previous fit results (=1) or not (=0);
seed     = 10;                                                             % Seed to generate random numbers;

if strcmp(noise, 'PN')
    itraj = [10 50 100];                                                   % Trajectories of the BD process to load;
elseif strcomp(noise, 'MNHo')
    sd = 0.5;                                                              % Standard deviation of the measurement error (log scale);
else
    var_a = 1;                                                             % parsameters of variance for MNHe;
    var_b = 2;
end

% ----------------------------------------------------------------------- %
% Values of model parsameters:
bS      = 0.63;                                                            % Birth rate of S in absence of antimicrobial;                         
bR      = 0.36;                                                            % Birth rate of R;

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

% ----------------------------------------------------------------------- %
% Time discretisation:
t0   = 0;                                                                  % Initial simulation time;
tf   = 48;                                                                 % Final simulation time;
ht   = 1e-3;                                                               % Time step for model simulation (solve ODEs);
texp = [2 4 6 8 10 12 16 20 24 36 48];                                     % Sampling times;

% ----------------------------------------------------------------------- %
% Discretisation of AMR level:
equi = 0;                                                                  % Equispaced discretisation of AMR level (=1) or not (=0)
                                                                           % If equi = 0 the user must define the AMR level discretisation
                                                                           % as a column array. For example: r = [0.3;0.35;0.7;0.99];
r    = [];                                                                 % Discretisation of AMR level if equi = 0;

if equi == 1 
    ra   = 0;                                                              % Minimum AMR level;
    rb   = 1;                                                              % Maximum AMR level;
    nr   = 50;                                                             % Number of subpopulations;
elseif numel(r) == 0
    fprintf('\n >> Please, define a discretisation of the AMR level or set "equi" option equal to one.\n')
    return
end

% Set ODE solver precision:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user-defined options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results_name = sprintf('Results/resPE_%s_%utraj.mat', noise, Ntraj);
    
rng(seed)                                                                  % Set seed for generate randoms;

pars = [bS;bR;alpha_b;d_maxS;alpha_d;beta_d;EC_50d;H_d;xi_SR;...            % Array with model parsameters;
       k_xi;N_T0;lambda_T0];

tmod     = t0:ht:tf;                                                       % Array with times to simulate model;
tmod     = sort(unique([tmod texp])).'; 
texp_ind = find(ismember(tmod, texp));  

if equi == 1                                                               % Array with AMR levels (if not previously defined);
    r  = linspace(ra, rb, nr)';
else
    nr = numel(r);
end

R = repmat(r, 1, nr) - repmat(r.', nr, 1);

nt    = numel(tmod);    
ntexp = numel(texp);
Nexp  = numel(Cexp);

% ----------------------------------------------------------------------- %
% Generate syntetic data for MNHo and MNHe if not previously generated
% (for PN case the data must be previously generated using SSA)

if load_res == 1
    
    load(results_name, 'NT_ave_data', 'pars_opt');
    
else
    
    % ------------------------------------------------------------------- %
    % Load data of SSA trajectories for the PN case:
    if strcmp(noise, 'PN')
        NT_data = zeros(ntexp, Nexp, Ntraj);
        for ii = itraj
            traj_name = sprintf('2024-PaperHeteroAMR/SSA/Results/resSSA_%03u.mat', itraj);
            load(traj_name, 'CFUST')
            NT_data(1:ntexp, 1:Nexp, ii) = CFUST(texp_ind, 1:Nexp);
        end

    else
        
        % --------------------------------------------------------------- %
        % Generate random trajectories for the MNHo or MNHe cases:
   
        f0 = exp(-lambda_T0*r);
        f0 = f0/sum(f0);
        N0 = N_T0*f0;                                                      % Initial condition for ODEs;
   
        NT = zeros(nt, Nexp);                                              % Initialice total population size;
                                                            
        R = R - triu(R) + tril(R).';
        Xi = xiSR*exp(kxi*(1 - R));
        Xi = Xi - diag(diag(Xi));

        AA_aux = Xi' - diag(sum(Xi, 2));                                   % Initialice coefficient matrix;  

        b     = bS*bR./(bR + r.^alpha_b*(bS - bR));                        % Birth rate;
        d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./...
            (beta_d^alpha_d + r.^alpha_d);                                 % Maximal kill rate;
     
        for iexp = 1:Nexp
    
            C  = Cexp(iexp);  
            HC = C^H_d/(C^H_d + EC_50d^H_d);
            d  = d_max*HC;                                                 % Death rate;

            AA = AA_aux + diag(b - d);                                     % Coefficient matrix for ODEs;

            [~, xout] = ode15s(@(t,s) Odes_cte(t, s, AA), tmod, N0, ODEoptions);
   
            
            NT(1:nt, iexp) = sum(xout, 2);                                 % Total population size;
        end

        % Add noise to data:
        NT_aux = NT(texp_ind, 1:Nexp);
        if strcmp(noise, 'MNHo')
            NT_data = log10(repmat(NT_aux, 1, 1, Ntraj)) + sd*randn(ntexp, Nexp, Ntraj);
        else
            var     = var_a*NT_aux.^var_b;
            NT_data = repmat(NT_aux, 1, 1, Ntraj) + repmat(sqrt(var), 1, 1, Ntraj).*randn(ntexp, Nexp, Ntraj);     
        end   
    end
    
    NT_ave_data = sum(NT_data, 3)/Ntraj;
    
end

% ----------------------------------------------------------------------- %
% Calibrate model with ESS:

bSmin        = 0.5;                                                        % Bounds for the (common) decision variables;
bSmax        = 5.0;

bRmin        = 1e-1;                                                         
bRmax        = bSmin - 1e-2;

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

if nr > 2
    problem.x_L = [bSmin;bRmin;alpha_bmin;d_maxSmin;alpha_dmin;beta_dmin;EC_50dmin;H_dmin;xi_SRmin;k_ximin;N_T0min;lambda_T0min]; 
    problem.x_U = [bSmax;bRmax;alpha_bmax;d_maxSmax;alpha_dmax;beta_dmax;EC_50dmax;H_dmax;xi_SRmax;k_ximax;N_T0max;lambda_T0max];
else
    problem.x_L = [bSmin;bRmin;d_maxSmin;EC_50dmin;H_dmin;xi_SRmin;N_T0min;lambda_T0min]; 
    problem.x_U = [bSmax;bRmax;d_maxSmax;EC_50dmax;H_dmax;xi_SRmax;N_T0max;lambda_T0max];
end

if strcmp(noise, 'MNHe')
    var_bmin    = 0;
    var_bmax    = 2;
    problem.x_L = [problem.x_L;var_bmin];
    problem.x_U = [problem.x_U;varb_max];
end
pars_nom    = problem.x_L + (problem.x_U - problem.x_L).*rand(size(problem.x_L));
problem.x_L = problem.x_L./pars_nom;
problem.x_U = problem.x_U./pars_nom;
np          = numel(problem.x_L);

if load_res == 1                                                           % Initial guess;
    problem.x_0 = pars_opt./pars_nom;
else
    problem.x_0 = problem.x_L + (problem.x_U - problem.x_L).*rand(size(problem.x_L));
end

problem.f = sprintf('costFun_%s',noise);                                   % Name of the cost function; 

opts.maxtime      = 1.0e6;                                                 % ESS solver options;
opts.maxeval      = 1.0e4;
opts.strategy     = 1;
opts.local.solver = 'fminsearch'; %'solnp';  %'wdn2fb';  %'fsqp';
opts.local.finish = 'fminsearch'; %'wdn2fb'; %'fsqp';
% opts.local.n1 = 50;
% opts.local.n2 = 50;

Results = ess_kernel(problem, opts, r, R, tmod, texp_ind, Cexp,...        % Call the optimization function (ESS):
                     NT_data, pars_nom, Weights, ODEoptions);

pars_opt = Results.xbest.*pars_nom;
f_best   = Results.fbest;
                
delete ess_report.mat                                                      % Remove the mat file generated by ESS

% ----------------------------------------------------------------------- %
% Save results:
if strcmp(noise, 'MNHo')
    save(results_name, 'r', 'tmod', 'texp', 'pars', 'pars_opt', 'f_best',...
         'NT_data','NT_ave_data','seed', 'sd')
elseif strcmp(noise, 'MNHe')
    save(results_name, 'r', 'tmod', 'texp', 'pars', 'pars_opt', 'f_best',...
         'NT_data','NT_ave_data','seed', 'var_a', 'var_b')
else
    save(results_name, 'r', 'tmod', 'texp', 'pars', 'pars_opt', 'f_best',...
         'NT_data','NT_ave_data','seed')
end

% ----------------------------------------------------------------------- %
% Plot results:

col1 = [39, 183, 222]/256;                                                  % Define colors,
col2 = [250, 128, 114]/256;
col3 = [75, 92, 56]/256;
col4 = [250, 141, 34]/256;
col5 = [90, 14, 45]/256;
col6 = [242, 220, 35]/256;

col  = [col1;col2;col3;col4;col5;col6];

mks{1} = 's';                                                              % Define markers,
mks{2} = 'diamond';
mks{3} = 'v';
mks{4} = '^';
mks{5} = '>';
mks{6} = '<';

if strcmp(noise, 'MNHo')
    logNT_ave_data = NT_ave_data;
else
    logNT_ave_data = log10(NT_ave_data);
end

PlotPE(tmod, r, par_opt, logNT_ave_data, col, mks)

rmpath('Functions')
rmpath('Results')



