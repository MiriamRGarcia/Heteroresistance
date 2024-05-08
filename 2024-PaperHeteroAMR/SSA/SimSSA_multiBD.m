%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSA_multiBD: SSA algorithm for the multivariate BD heteroresistance model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CFUS, CFUST] = SimSSA_multiBD(tmod, r, pars, Cexp, N_TL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% tmod  = Array of simulation times to almacenate BD trajectories (nt x 1);
% r     = Discretisation of the AMR level (values between entire sensitivity
%         S = 0 and entire resistance R = 1) (nr x 1),
% pars  = Array of model parameters,
% Cexp  = Antimicrobial concentrations (constant),
% N_TL  = Maximum cell count to stop simulation,
% seed  = Seed to generate random numbers for reproducibility,
%
% OUTPUT:
% CFUS  = Array with the multivariate BD process at the model times tmod
%         (size nt x nr x Nexp),
% CFUST = Array with the total counts (nt + Nexp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nexp = numel(Cexp);

% ----------------------------------------------------------------------- %
% Calculate BDM rates:

% Obtain parameter values:
bS        = pars(1);                                                      
bR        = pars(2);

alpha_b   = pars(3);

d_maxS    = pars(4);
alpha_d   = pars(5);
beta_d    = pars(6);

EC_50d    = pars(7);
H_d       = pars(8);

xi_SR     = pars(9);
k_xi      = pars(10);

N_T0      = pars(11);
lambda_T0 = pars(12);

% Matrix with modification rates:
nr  = numel(r);                                                         % Size of the AMR level discretisation,   

RR  = repmat(r, 1, nr) - repmat(r.', nr, 1);                         % Auxiliary matrix to calculate modification rate,
RR  = RR - triu(RR) + tril(RR).';

Xi  = xi_SR*exp(k_xi*(1 - RR));                                            % Matrix with modification rates,                                       
Xi  = Xi - diag(diag(Xi));

% Birth rate:
b    = bR*bR./(bR + r.^alpha_b*(bS - bR));            

% Maximal death rate:
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);                                                    


% ----------------------------------------------------------------------- %
% Preeliminary calculations for SSA:

% Time discretisation:
nt = numel(tmod);                                                          % Size of the time discretisation,
t0 = tmod(1);                                                              % Initial simulation time,
tf = tmod(nt);                                                             % Final simulation time,

% Initial heteroresistance distribution:
f0 = exp(-lambda_T0*r);                                                    % Exponentially decaying subpopulations frequencies,
f0 = f0/sum(f0);

N0    = floor(N_T0*f0);                                                    % Initial number of cells,
N0(1) = N0(1) + N_T0 - sum(N0);

% Calculate matrix of state transitions:
trans    = [];
can_bas  = eye(nr, nr);

for jj = 1:nr
    trans_aux = zeros(nr, nr);
    
    for ii = 1:nr
        trans_aux(ii, 1:nr) = can_bas(1:nr, jj).' - can_bas(1:nr, ii).';
    end
    trans_aux(jj, 1:nr) = can_bas(1:nr, jj).';
    trans_aux = [trans_aux;-can_bas(1:nr, jj).'];
    trans     = [trans;trans_aux];
end

% Number of possible transitions:
n_re = size(trans, 1);

% Auxiliary matrix for propensities of BM reactions:
prp_BM = diag(b) + Xi;

% ----------------------------------------------------------------------- %
% Run SSA:
CFUS  = zeros(nt, nr, Nexp);    
   
for iexp = 1:Nexp
    
    % Initial condition:
    CFUS(1, 1:nr, iexp) = N0;
    
    % Kill rate:
    HC = Cexp(iexp)^H_d/(Cexp(iexp)^H_d + EC_50d^H_d);
    d  = d_max*HC;  

    % ------------------------------------------------------------------- %
    % Initialise variables:
    tsim        = t0;                                                      % Initialise simulation time,
    tint_count  = 1;                                                       % Initialise counter for the current interval of time discretisation,

    cell_counts = N0;                                                      % Initialise cell counts,

    while tsim < tf                                   
    
        % --------------------------------------------------------------- %
        % Calculate propensity function:
        prp_b = prp_BM.*repmat(cell_counts, 1, nr);
        prp_d = (d.*cell_counts).';
    
        prp   = [prp_b;                                                        % Array of propensities for reactions,
                 prp_d];
         
        prp   = reshape(prp, [], 1);

        % --------------------------------------------------------------- %
        % Direct SSA method:
        prp_sum = sum(prp);                                                % Sum of propensities,
    
        tau  = -log(rand)/prp_sum;                                         % Calculate time to the next reaction,
        tsim = tsim + tau;                                                 % Actualice simulation time,

    
        aux = rand*prp_sum;                                                % Calculate index of the next reaction,
        isum = 0;
        for ii = 1:n_re
            isum = isum + prp(ii);
            if isum > aux
                inext = ii;
                break
            end
        end
    
        %---------------------------------------------------------------- %
        % Almacenate the BD trajectories in the discretisation times:
        tind = find(tmod > tsim);
    
        if tind - 1 > tint_count                                           % Comprobate if discretisation interval changed,
            tint_count = tint_count + 1;
            CFUS(tint_count, 1:nr, iexp) = cell_counts;
        end  
    
        cell_counts = cell_counts + trans(inext, 1:nr).';
    
        N_T = sum(cell_counts);                                            % Total cell counts,

        if N_T > N_TL                                                      % Stop simulation if total counts reach the limit value,
        
            tint_count = tint_count + 1;
            CFUS(tint_count, 1:nr, iexp) = cell_counts;
            CFUS(tint_count + 1:nt, 1:nr, iexp) = NaN;
            break
        
        elseif N_T < 1                                                     % Stop simulation if the population goes extinct,
        
            CFUS(tint_count+1:nt, 1:nr, iexp) = zeros(nt - tint_count, 1:nr);
            break
        
        end

    end

    % ----------------------------------------------------------------------- %
    % If the time to next reaction exceds the final time:
    CFUS(tint_count:nt, 1:nr, iexp) = repmat(CFUS(tint_count, 1:nr, iexp), nt - tint_count + 1, 1);

end

CFUST = sum(CFUS, 2);

end