%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SimSSA_multiBD: SSA (direct method) to generate trajectories of
%                 the multivariate BD heteroresistance model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CFUS, CFUST] = SimSSA_multiBD(tmod, r, pars, Cexp, N_TL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% tmod  = Array of simulation times to almacenate BD trajectories (nt x 1);
% r     = Discretisation of the AMR level (values between entire sensitivity
%         S = 0 and entire resistance R = 1) (nr x 1),
% pars  = Array of model parameters,
% Cexp  = Antimicrobial concentrations (assumed constant in each experiment),
% N_TL  = Maximum cell count to stop simulation,
% seed  = Seed to generate uniform random numbers (for data reproducibility),
%
% OUTPUT:
% N   = Array with the multivariate BD process at the model times tmod
%         (size: nt x nr x Nexp),
% N_T = Array with the total counts (size: nt x Nexp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nt   = numel(tmod);
nr   = numel(r);                                                    
Nexp = numel(Cexp);

% ----------------------------------------------------------------------- %
% Calculate BDM rates:

% Obtain parameter values:
b_S       = pars(1);                                                      
b_R       = pars(2);

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

% Build matrix with modification rates:
RR  = repmat(r, 1, nr) - repmat(r.', nr, 1);                   
RR  = RR - triu(RR) + tril(RR).';

Xi  = xi_SR*exp(k_xi*(1 - RR));                                                                        
Xi  = Xi - diag(diag(Xi));

% Birth rate:
b     = b_S*b_R./(b_R + r.^alpha_b*(b_S - b_R));            

% Maximal death rate:
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);                                                    

% ----------------------------------------------------------------------- %
% Preeliminary calculations for SSA:

% Time discretisation:
t0 = tmod(1);                                                             
tf = tmod(nt);                                                        

% Initial heteroresistance distribution:
f0 = exp(-lambda_T0*r);                                                 
f0 = f0/sum(f0);

N_0    = floor(N_T0*f0);                                               
N_0(1) = N0(1) + N_T0 - sum(N0);

% Calculate matrix with state transitions:
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
N  = zeros(nt, nr, Nexp);    
   
for iexp = 1:Nexp
    
    % Initial condition:
    N(1, 1:nr, iexp) = N_0;
    
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
        % Calculate reaction rates:
        prp_b = prp_BM.*repmat(cell_counts, 1, nr);
        prp_d = (d.*cell_counts).';
    
        prp   = [prp_b;                                                    % Array of reaction rates,
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
            N(tint_count, 1:nr, iexp) = cell_counts;
        end  
    
        cell_counts = cell_counts + trans(inext, 1:nr).';
    
        N_T = sum(cell_counts);                                            % Total cell counts,

        if N_T > N_TL                                                      % Stop simulation if total counts reached the limit value,
        
            tint_count = tint_count + 1;
            N(tint_count, 1:nr, iexp) = cell_counts;
            N(tint_count + 1:nt, 1:nr, iexp) = NaN;
            break
        
        elseif N_T < 1                                                     % Stop simulation if the population becomes extinct,
        
            N(tint_count+1:nt, 1:nr, iexp) = zeros(nt - tint_count, 1:nr);
            break
        
        end

    end

    % ----------------------------------------------------------------------- %
    % If the time to next reaction exceds the final time:
    N(tint_count:nt, 1:nr, iexp) = repmat(N(tint_count, 1:nr, iexp), nt - tint_count + 1, 1);

end

% Total cell counts:
N_T = sum(N, 2);

end
