%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SensMultiExp: Sensitivities of the average heteroresistance model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N_T, sensN_T] = SensMultiExp(tmod, r, C, pars, N_0, ODEoptions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% tmod       = Simulation times (nt x 1);
% r          = Discretisation of AMR level (nr x 1);
% C          = Array of (constant) antimicrobial concentrations (Nexp x 1);
% pars       = Values of model parameters (np x 1);
% N_0        = Initial cell counts (nr x 1);
% ODEoptions = Options to solve ODEs;
%
% OUTPUT:
% N_T     = Average cell counts (nt x Nexp),
% sensN_T = Sensitivity of cell counts to model parameters (nt x Nexp x np);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% Problem sizes:
nC = numel(C);                                                     
nt = numel(tmod);                                                    
nr = numel(r);                                                        
np = numel(pars);                                                   

% ----------------------------------------------------------------------- %
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

% ----------------------------------------------------------------------- %
% Initial condition of the sensitivity system (same for each experiment):

% Initial condition for average and sensitivities:
s0 = [N_0;
      zeros(nr*(np - 2),1)];                                               % Initial conditions are independent of parameters
                                                                           % unless for X0 and lamb_IC
% Ampliate with the initial condition for sensitivity to IC parameters:
f0_aux     = exp(-lambda_T0*r);
s0_N0      = f0_aux/sum(f0_aux);
RR_aux     = repmat(r.', nr, 1) - repmat(r, 1, nr);    
s0_lamb_IC = N_T0/(sum(f0_aux)^2)*f0_aux.*(RR_aux*f0_aux);


s0 = [s0;
      s0_N0;
      s0_lamb_IC];

% ----------------------------------------------------------------------- %
% Common calculations for ODEs (independent of the experiment):

% Calculate mutation-rates matrix:
RR  = repmat(r, 1, nr) - repmat(r.', nr, 1);                         
RR  = RR - triu(RR) + tril(RR).';

Xi = xi_SR*exp(k_xi*(1 - RR));
Xi = Xi - diag(diag(Xi));

% Auxiliary coefficient matrix:
AA_aux = Xi - diag(sum(Xi, 2));

% Derivative of the dynamics with respect to bS:
db_S(1)  = 1;
db_S(nr) = 0;
db_S(2:(nr-1)) = b_R^2*(1 - r(2:(nr-1)).^alpha_b)./(r(2:(nr-1)).^alpha_b*(b_R - b_S) - b_R).^2;

AA_bS = diag(db_S);

% Derivative with respect to bR:
db_R(1)  = 0;
db_R(nr) = 1;
db_R(2:(nr-1)) = b_S^2*r(2:(nr-1)).^alpha_b./(b_R + r(2:(nr-1)).^alpha_b*(b_S - b_R)).^2;

AA_bR    = diag(db_R);

% Derivative with respect to alpha_b:
dalpha_b(1)  = 0;
dalpha_b(nr) = 0;
dalpha_b(2:(nr-1)) = b_R*b_S*(b_R - b_S)*log(r(2:(nr-1))).*r(2:(nr-1)).^alpha_b./(b_R + r(2:(nr-1)).^alpha_b*(b_S - b_R)).^2;

AA_alpha_b  = diag(dalpha_b);

% Derivative with respect to xi_SR:
dxi_SR  = exp(k_xi*(1 - RR));
dxi_SR  = dxi_SR - diag(diag(dxi_SR));
AA_xiSR = dxi_SR - diag(sum(Xi/xi_SR, 2));

% Derivative with respect to k_xi:
dk_xi  = (1 - RR).*Xi;
dk_xi  = dk_xi - diag(diag(dk_xi));
AA_kxi = dk_xi - diag(sum(dk_xi, 2));

% Calculate birth rate:
b     = b_S*b_R./(b_R + (b_S - b_R)*r.^alpha_b);
 
% Calculate maximal kill rate:
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);

% Solve ODEs sensitivity system for each experiment:
N_T         = zeros(nt, nC);
sensN_T    = zeros(nt, np, nC);

dalpha_d  = zeros(nr,1);
dbeta_d   = zeros(nr,1);
dH_d      = zeros(nr,1);
dd_maxS   = zeros(nr,1);
dEC_50d   = zeros(nr,1);

for iexp = 1:nC
    
    %-------------------------------------------------%
    % Calculate coefficient matrix of the state system:
    
    % Drug concentration:
    CC = C(iexp);
    
    % Calculate kill rate at current time:
    HC = CC^H_d/(CC^H_d + EC_50d^H_d);
    d  = d_max*HC;

    AA = AA_aux + diag(b - d);
    
    %-------------------------------------------------%
    % Calculate derivatives of the coefficient matrix:
    
     % Derivative of A with respect to d_maxS:
    dd_maxS(1)  = 1;
    dd_maxS(nr) = 0;
    dd_maxS(2:(nr - 1)) = beta_d^alpha_d*(1 - r(2:(nr - 1)).^alpha_d)./(beta_d^alpha_d + r(2:(nr - 1)).^alpha_d);

    dd_maxS  = - HC*dd_maxS;

    AA_dmaxS = diag(dd_maxS);
    
    % Derivative with respect to alph_d:
    dalpha_d(1)  = 0;
    dalpha_d(nr) = 0;
    dalpha_d(2:(nr -1)) = beta_d^alpha_d*d_maxS*r(2:(nr -1)).^alpha_d.*((1 + beta_d^alpha_d)*log(r(2:(nr - 1))) + (r(2:(nr -1)).^alpha_d - 1)*log(beta_d))./(r(2:(nr -1)).^alpha_d + beta_d^alpha_d).^2;

    dalpha_d  = HC*dalpha_d;

    AA_alpha_d = diag(dalpha_d);
    
    % Derivative with respect to beta_d:
    dbeta_d(1)  = 0;
    dbeta_d(nr) = 0;
    dbeta_d(2:(nr -1)) =  alpha_d*beta_d^(alpha_d - 1)*d_maxS*(r(2:(nr -1)).^alpha_d - 1).*r(2:(nr -1)).^alpha_d./(beta_d^alpha_d + r(2:(nr -1)).^alpha_d).^2;

    dbeta_d  = HC*dbeta_d;
    
    AA_beta_d = diag(dbeta_d);

    % Derivative with respect to EC50k:
    dEC_50d(nr)         = 0;
    dEC_50d(1:(nr - 1)) = - EC_50d^(H_d - 1)*H_d*beta_d^alpha_d*d_maxS*(r(1:(nr -1)).^alpha_d - 1)./(r(1:(nr -1)).^alpha_d + beta_d^alpha_d);
    dEC_50d             = HC^2*dEC_50d;

    AA_EC_50d            = diag(dEC_50d);

    % Derivative with respect to Hk:
    if CC > 0
        dH_d(1:(nr -1)) = EC_50d^H_d*beta_d^alpha_d*d_maxS*(log(CC) - log(EC_50d))*(r(1:(nr -1)).^alpha_d - 1)./(r(1:(nr -1)).^alpha_d + beta_d^alpha_d);
        dH_d(nr)        = 0;
        dH_d            = HC^2*dH_d;
    else
        dH_d(1:nr)      = 0;
    end
    
    AA_H_d               = diag(dH_d);
    
    %-------------------------------------------------%
    % Calculate output sensitivities:
    [~, St_sens] = ode15s(@(t,s) OdesSens(t, s, AA, AA_bS, AA_bR, AA_alpha_b,...
                   AA_dmaxS, AA_alpha_d, AA_beta_d, AA_EC_50d, AA_H_d, AA_xiSR, AA_kxi), tmod, s0, ODEoptions);
   
    % Total population:
    N              = St_sens(1:nt, 1:nr);
    N_T(1:nt, iexp) = sum(N, 2);
      
    % Sensitivities of the states:
    sensN_bS        = St_sens(1:nt,nr + 1:2*nr);
    sensN_bR        = St_sens(1:nt,2*nr + 1:3*nr);
    sensN_alpha_b   = St_sens(1:nt,3*nr + 1:4*nr);
    sensN_d_maxS    = St_sens(1:nt,4*nr + 1:5*nr);
    sensN_alpha_d   = St_sens(1:nt,5*nr + 1:6*nr);
    sensN_beta_d    = St_sens(1:nt,6*nr + 1:7*nr);
    sensN_EC_50d    = St_sens(1:nt,7*nr + 1:8*nr);
    sensN_H_d       = St_sens(1:nt,8*nr + 1:9*nr);
    sensN_xiSR      = St_sens(1:nt,9*nr + 1:10*nr);
    sensN_kxi       = St_sens(1:nt,10*nr + 1:11*nr);
    sensN_NT0       = St_sens(1:nt,11*nr + 1:12*nr);
    sensN_lambda_T0 = St_sens(1:nt,12*nr + 1:13*nr);
    
    % Sensitivity of the total population:
    sensN_T(1:nt, 1, iexp)  = sum(sensN_bS, 2);
    sensN_T(1:nt, 2, iexp)  = sum(sensN_bR, 2);
    sensN_T(1:nt, 3, iexp)  = sum(sensN_alpha_b, 2);
    sensN_T(1:nt, 4, iexp)  = sum(sensN_d_maxS, 2);
    sensN_T(1:nt, 5, iexp)  = sum(sensN_alpha_d, 2);
    sensN_T(1:nt, 6, iexp)  = sum(sensN_beta_d, 2);
    sensN_T(1:nt, 7, iexp)  = sum(sensN_EC_50d, 2);
    sensN_T(1:nt, 8, iexp)  = sum(sensN_H_d, 2);
    sensN_T(1:nt, 9, iexp)  = sum(sensN_xiSR, 2);
    sensN_T(1:nt, 10, iexp) = sum(sensN_kxi, 2);
    sensN_T(1:nt, 11, iexp) = sum(sensN_NT0, 2);
    sensN_T(1:nt, 12, iexp) = sum(sensN_lambda_T0, 2);
    
end

          
end