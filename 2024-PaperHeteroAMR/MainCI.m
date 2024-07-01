%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MainCI: Main file to calculate FIM-based confidence intervals 
%         for calibration of the BD heteroresistance model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBSERVATION:
% To calculate the FIM-based confidence intervals for parameter estimates
% obtained through Maximum Likelihood Estimation, previous 
% calibration results (generated with MainPE.m) must be available 
% at folder: Results/ResPE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables

% Add folder with necessary functions:
addpath('Functions')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Confidence level:
confLev    = 0.9;

% Noise assumption (= 'MNHo'; 'MNHe'; 'PN');
noise      = 'MNHe';   
    
% Name of the file keeping the calibration results:
load_name  = 'ResPE_50subpop_MNHe_3traj_run1.mat';

% Set ODEs solver (ode15s) precision to calculate state sensitivities:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user-defined settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
fprintf('\n>> The user has selected to load calibration results from: %s', load_name)
fprintf('\n>> The user has selected the noise assumption: %s', noise)

if numel(strfind(load_name, noise)) < 1
    fprintf('\n>> The selected noise assumption does not coincide with that used for calibration.')
    fprintf('\n>> Please, change the name of the file with calibration results or the noise assumption and run again.\n')
    return
end
full_load_name  = sprintf('Results/ResPE/%s', load_name);

% Initialice auxiliary array with parameters of variance in MNHe case:
pars_var  = [];

if exist(full_load_name, 'file')
    
    % Load optimal parameters and data from previous calibration results:
    if strcmp(noise, 'MNHo')
        load(full_load_name, 'r', 'tsim', 'Cexp', 'texp', 'pars_opt',...
            'sd_opt', 'Weights', 'Var_data')
        
    elseif strcmp(noise, 'MNHe')
        load(full_load_name, 'r', 'tsim', 'Cexp', 'texp', 'pars_opt',...
            'Weights', 'Var_data', 'var_a_opt', 'var_b_opt')
        
        % Remove optimal var_b of the model parameters:
        pars_opt = pars_opt(1:end - 1);
        
        % Optimal parameters for variance:
        pars_var = [var_a_opt;var_b_opt];
        
    else
        load(full_load_name, 'r', 'tsim', 'Cexp', 'texp', 'pars_opt',...
            'Weights', 'Var_data')
    end

else
    
    fprintf('\n>> The file: %s cannot be found in folder: Results/resPE', load_name)
    fprintf('\n>> This script, however, needs previous calibration results to work.')
    fprintf('\n>> Please, check the name of the file defined in: load_name, and run again.')
    
    return
end


% Indexes of sampling times within simulation times:
texp_ind = find(ismember(tsim, texp));

% Problem sizes:
m_r    = numel(r);
m_t    = numel(tsim);    
m_texp = numel(texp);
m_p    = numel(pars_opt);
m_e    = numel(Cexp);

% ----------------------------------------------------------------------- %
% Calculate sensitivity matrix:

fprintf('\n>> Calculating sensitivities...')

% Initial condition:
f0  = exp(- pars_opt(end)*r);
f0  = f0/sum(f0);
N_0 = pars_opt(end - 1)*f0;

% Call function to calculate state sensitivities of the average BD model
% with the calibrated parameters:
[N_T, sensN_T] = SensME(tsim, r, Cexp, pars_opt, N_0, ODEoptions);

% Obtain average and sensitivities at the sampling times:
N_T     = N_T(texp_ind, 1:m_e);
sensN_T = sensN_T(texp_ind, 1:m_p, 1:m_e);

% ----------------------------------------------------------------------- %
% Build covariance matrix:
if strcmp(noise, 'MNHo')
    
    CovMatrix = sd_opt^2*ones(m_texp*m_e, m_texp*m_e);
    
elseif strcmp(noise, 'MNHe')
    
    auxN_T    = reshape(N_T.', [], 1); 
    CovMatrix = diag(var_a_opt*auxN_T.^var_b_opt);
    
else
    
    % Remove NaN data:
    NaN_ind           = find(Weights == 0);
    Var_data(NaN_ind) = NaN;
    
    Var_data          = reshape(Var_data.', [], 1);
    
    Var_data(isnan(Var_data)) = [];
    
    % Build covariance matrix with sample variances:
    CovMatrix = diag(Var_data);
    
end

% ----------------------------------------------------------------------- %
% Build sensitivity matrix:
fprintf('\n>> Building sensitivity matrix...')

SensMatrix = [];

if strcmp(noise, 'MNHo')
    
    % Log10 of total average counts:
    y_T = log10(N_T);

    % Calculate sensitivities of y_T to parameters:
    sensy_T = zeros(m_texp, m_p, m_e);

    for ip = 1:m_p

        aux = reshape(sensN_T(1:m_texp, ip, 1:m_e), m_texp, m_e);
        
        % Obtain sensitivities at the sampling times (log scale):
        aux = (aux./N_T)/log(10);
        
        % Almacenate sensitivity of y_T = log(N_T):
        sensy_T(1:m_texp, ip, 1:m_e) = aux;
        
        % Log-scaling of parameters:
        aux = pars_opt(ip)*aux;
        
        SensMatrix = [SensMatrix reshape(aux.', [], 1)];
    end
else
    for ip = 1:m_p
    
        % Obtain sensitivities at the sampling times:
        aux = reshape(sensN_T(1:m_texp, ip, 1:m_e), m_texp, m_e);

        % Log scaling of parameters:
        aux = pars_opt(ip)*aux;
        
        if strcmp(noise, 'PN')
            aux(NaN_ind) = NaN;
        end
        
        aux             = reshape(aux.', [], 1);
        
        aux(isnan(aux)) = [];
        
        SensMatrix = [SensMatrix aux];
    end
end

% Normalised sensitivity matrix using the sampling times:      
% Important! this is the normalised sensitivity with respect parameters
% before log-scaling of parameters! 
normSensMatrix = [];

if strcmp(noise, 'MNHo')
     for ip = 1:m_p
        aux = reshape(sensy_T(1:m_texp, ip, 1:m_e), m_texp, m_e);
        aux = pars_opt(ip)*aux./y_T;
        normSensMatrix = [normSensMatrix reshape(aux.',[],1)];
    end   
else
    for ip = 1:m_p
        aux = reshape(sensN_T(1:m_texp, ip, 1:m_e), m_texp, m_e);
        aux = pars_opt(ip)*aux./N_T;
        
        if strcmp(noise, 'PN')
            aux(NaN_ind) = NaN;
        end
        
        aux             = reshape(aux.', [], 1);
        
        aux(isnan(aux)) = [];
        normSensMatrix  = [normSensMatrix aux];
    end
end

%%
% ----------------------------------------------------------------------- %
% Calculate confidence intervals from FIM:

fprintf('\n>> Calculating the FIM and confidence intervals...')

[FIM, CI] = FIM_CI(pars_opt, pars_var, SensMatrix, normSensMatrix,...
                   CovMatrix, confLev, noise, m_r, m_texp, m_e); 

% ----------------------------------------------------------------------- %
% Save results:
res_name = sprintf('Results/ResFIM/FIM_%s', load_name);
save(res_name, 'r', 'tsim', 'texp', 'Cexp', 'pars_opt', 'pars_var', 'CovMatrix',...
    'SensMatrix', 'normSensMatrix', 'FIM', 'CI')

% Remove folder with functions from path:
rmpath('Functions')

