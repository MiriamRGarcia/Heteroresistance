%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MainSSA_multiBD: Main file to simulate trajectories of the 
%                  multivariate BD heteroresistance model with
%                  constant antimicrobial concentration using the SSA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined settings:

% Choose implementation (direct method = SSA or rejection based = RSSA):
method = 'RSSA'; % = 'SSA'; = 'RSSA';

% Number of trajectories:
m_traj = 10;

% Threshold on total cell counts to stop simulation if neccesary:
N_TL   = 1e7;

% Define parameter values:
b_S       = 0.63;                                                          % Natural birth rate of S cells;
b_R       = 0.36;                                                          % Natural birth rate of R cells;
alpha_b   = 2;                                                             % Shape coefficient of AMR fitness cost;

d_maxS    = 3.78;                                                          % Maximal kill rate of S cells;
alpha_d   = 3;                                                             % Shape coefficient of death rate with AMR level;
beta_d    = 0.4;                                                           % Shape coefficient of death rate with AMR level;

EC_50d    = 1;                                                             % Half maximal inhibitory concentration;
H_d       = 1;                                                             % Hill coefficient of drug inhibition;

xi_SR     = 1e-6;                                                          % Modification rate in AMR level between S and R cells;
k_xi      = log(1e2);                                                      % Decay velocity of modification rate with AMR level;

N_T0      = 1e6;                                                           % Initial total cell count;
lambda_T0 = 50;                                                            % Decay velocity of initial heteroresistance distribution with AMR level;

pars      = [b_S;b_R;alpha_b;d_maxS;alpha_d;beta_d;EC_50d;H_d;...          % Parameter array;
             xi_SR;k_xi;N_T0;lambda_T0];


% Array of antimicrobial concentration (assumed constant in each experiment):
MIC_S = EC_50d*(b_S/(d_maxS - b_S))^(1/H_d);                               % Minimum inhibitory concentration of S cells,
Cexp  = MIC_S*[0 1 2 4 8 32].';                                            % Array of antimicrobial concentration,

% Time discretisation:
t0   = 0;                                                                  % Initial time [h];
tf   = 48;                                                                 % Final time [h]; 
ht   = 1e-3;                                                               % Time step;        
tsim = t0:ht:tf;                                                           % Time discretisation;    

% Discretisation of AMR level:
ra   = 0;                                                                  % Minimum AMR level (between entire sensitivity = 0 and entire resistance = 1),
rb   = 1;                                                                  % Maximum AMR level (between entire sensitivity = 0 and entire resistance = 1),
m_r  = 50;                                                                 % Size of AMR level discretisation,
r    = linspace(ra, rb, m_r).';                                            % AMR level discretisation,

% End of user-defined settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% Preeliminary calculations:

% Problem sizes:
m_t       = numel(tsim);
m_e       = numel(Cexp);

% Name of the function implementing the method:
fun_name  = sprintf('multiBD_%s', method);

% Calculate matrix of state transitions:
can_basis = eye(m_r, m_r);
trans     = [can_basis -can_basis];                                        % Define possible state transitions;
 
for ii = 1:m_r
    trans_aux = zeros(m_r, m_r);

    for jj = 1:m_r
        trans_aux(1:m_r, jj) = can_basis(1:m_r, jj) - can_basis(ii, 1:m_r).';
    end
    
    trans_aux(:, ii) = [];
    
    trans = [trans trans_aux];
end

% Calculate BD rates in each experiment:
b     = b_S*b_R./(b_R + r.^alpha_b*(b_S - b_R));                           % Array of birth rates;
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);  %Array with maximal death rates;
d     = Cexp.^H_d./(Cexp.^H_d + EC_50d^H_d)*d_max.';                       % Matrix m_e x m_r with deaths rates in each experiment;

% Calculate modification rates in AMR level:
RR  = repmat(r, 1, m_r) - repmat(r.', m_r, 1);                           
RR  = RR - triu(RR) + tril(RR).';

Xi  = xi_SR*exp(k_xi*(1 - RR));                                                                        
Xi  = Xi - diag(diag(Xi));

% Calculate initial cell counts:
f_0    = exp(-lambda_T0*r);
f_0    = f_0/sum(f_0);
N_0    = floor(N_T0*f_0);
N_0(1) = N_0(1) + (N_T0 - sum(N_0));

% ----------------------------------------------------------------------- %
% Call to SSA and save results:
N   = zeros(m_t, m_r, m_e);                                                % Initialise trajectories of cell counts;
N_T = zeros(m_t, m_e);                                                     % Initialise trajectories of total cell count;

for itraj = 1:m_traj

    % Save seed for data reproducibility:
    seed  = itraj*100;
    
    % Name of the file to keep the results:
    res_name = sprintf('Results/res%s_%03u', method, itraj);
    
    for iexp = 1:m_e
        
        % Set same seed for the different experiments:
        rng(seed)
        
        % Call to tRSSA:
        [N_aux, N_T_aux] = feval(fun_name, tsim, b, d(iexp, 1:m_r).',...
                                 Xi, trans, N_0, N_TL);
    
        N(1:m_t, 1:m_r, iexp) = N_aux;
        N_T(1:m_t, iexp)      = N_T_aux;
    end
    
    save(res_name, 'r', 'tsim', 'pars', 'Cexp', 'seed', 'N_TL', 'N', 'N_T')
end