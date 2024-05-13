%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate sensitivities of the average heteroresistance model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NT, sens_NT] = SensMultiExp(tmod, r, C, par, N0, ODEoptions)

% ----------------------------------------------------------------------- %
% Problem sizes:
nC = numel(C);                                                     
nt = numel(tmod);                                                    
nr = numel(r);                                                        
np = numel(par);                                                   

% ----------------------------------------------------------------------- %
% Obtain parameter values:
bS        = par(1);
bR        = par(2);
alpha_b   = par(3);

d_maxS    = par(4);
alpha_d   = par(5);
beta_d    = par(6);

EC_50d    = par(7);
H_d       = par(8);

xi_SR     = par(9);
k_xi      = par(10);

N_T0      = par(11);
lambda_T0 = par(12);

% ----------------------------------------------------------------------- %
% Initial condition of the sensitivity system (same for each experiment):

% Initial condition for average and sensitivities:
s0 = [N0;
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
d_mugS(1)  = 1;
d_mugS(nr) = 0;
d_mugS(2:(nr-1)) = bR^2*(1 - r(2:(nr-1)).^alpha_b)./(r(2:(nr-1)).^alpha_b*(bR - bS) - bR).^2;

AA_bS = diag(d_mugS);

% Derivative with respect to bR:
d_mugR(1)  = 0;
d_mugR(nr) = 1;
d_mugR(2:(nr-1)) = bS^2*r(2:(nr-1)).^alpha_b./(bR + r(2:(nr-1)).^alpha_b*(bS - bR)).^2;

AA_bR    = diag(d_mugR);

% Derivative with respect to alpha_b:
d_alphg(1)  = 0;
d_alphg(nr) = 0;
d_alphg(2:(nr-1)) = bR*bS*(bR - bS)*log(r(2:(nr-1))).*r(2:(nr-1)).^alpha_b./(bR + r(2:(nr-1)).^alpha_b*(bS - bR)).^2;

AA_alpha_b  = diag(d_alphg);

% Derivative with respect to xi_SR:
d_xiSR  = exp(k_xi*(1 - RR));
d_xiSR  = d_xiSR - diag(diag(d_xiSR));
AA_xiSR = d_xiSR - diag(sum(Xi/xi_SR, 2));

% Derivative with respect to k_xi:
d_kxi  = (1 - RR).*Xi;
d_kxi  = d_kxi - diag(diag(d_kxi));
AA_kxi = d_kxi - diag(sum(d_kxi, 2));

% Calculate birth rate:
b     = bS*bR./(bR + (bS - bR)*r.^alpha_b);
 
% Calculate maximal kill rate:
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);

% Solve ODEs sensitivity system for each experiment:
NT         = zeros(nt, nC);
sens_NT    = zeros(nt, np, nC);

d_alpha_d  = zeros(nr,1);
d_beta_d   = zeros(nr,1);
d_H_d      = zeros(nr,1);
d_d_maxS   = zeros(nr,1);
d_EC_50d   = zeros(nr,1);

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
    d_d_maxS(1)  = 1;
    d_d_maxS(nr) = 0;
    d_d_maxS(2:(nr - 1)) = beta_d^alpha_d*(1 - r(2:(nr - 1)).^alpha_d)./(beta_d^alpha_d + r(2:(nr - 1)).^alpha_d);

    d_d_maxS  = - HC*d_d_maxS;

    AA_dmaxS = diag(d_d_maxS);
    
    % Derivative with respect to alph_d:
    d_alpha_d(1)  = 0;
    d_alpha_d(nr) = 0;
    d_alpha_d(2:(nr -1)) = beta_d^alpha_d*d_maxS*r(2:(nr -1)).^alpha_d.*((1 + beta_d^alpha_d)*log(r(2:(nr - 1))) + (r(2:(nr -1)).^alpha_d - 1)*log(beta_d))./(r(2:(nr -1)).^alpha_d + beta_d^alpha_d).^2;

    d_alpha_d  = HC*d_alpha_d;

    AA_alpha_d = diag(d_alpha_d);
    
    % Derivative with respect to beta_d:
    d_beta_d(1)  = 0;
    d_beta_d(nr) = 0;
    d_beta_d(2:(nr -1)) =  alpha_d*beta_d^(alpha_d - 1)*d_maxS*(r(2:(nr -1)).^alpha_d - 1).*r(2:(nr -1)).^alpha_d./(beta_d^alpha_d + r(2:(nr -1)).^alpha_d).^2;

    d_beta_d  = HC*d_beta_d;
    
    AA_beta_d = diag(d_beta_d);

    % Derivative with respect to EC50k:
    d_EC_50d(nr)         = 0;
    d_EC_50d(1:(nr - 1)) = - EC_50d^(H_d - 1)*H_d*beta_d^alpha_d*d_maxS*(r(1:(nr -1)).^alpha_d - 1)./(r(1:(nr -1)).^alpha_d + beta_d^alpha_d);
    d_EC_50d             = HC^2*d_EC_50d;

    AA_EC_50d            = diag(d_EC_50d);

    % Derivative with respect to Hk:
    if CC > 0
        d_H_d(1:(nr -1)) = EC_50d^H_d*beta_d^alpha_d*d_maxS*(log(CC) - log(EC_50d))*(r(1:(nr -1)).^alpha_d - 1)./(r(1:(nr -1)).^alpha_d + beta_d^alpha_d);
        d_H_d(nr)        = 0;
        d_H_d            = HC^2*d_H_d;
    else
        d_H_d(1:nr)      = 0;
    end
    
    AA_H_d               = diag(d_H_d);
    
    %-------------------------------------------------%
    % Calculate output sensitivities:
    [~, St_sens] = ode15s(@(t,s) OdesSens(t, s, AA, AA_bS, AA_bR, AA_alpha_b,...
                   AA_dmaxS, AA_alpha_d, AA_beta_d, AA_EC_50d, AA_H_d, AA_xiSR, AA_kxi), tmod, s0, ODEoptions);
   
    % Total population:
    N              = St_sens(1:nt, 1:nr);
    NT(1:nt, iexp) = sum(N, 2);
      
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
    sens_NT(1:nt, 1, iexp)  = sum(sensN_bS, 2);
    sens_NT(1:nt, 2, iexp)  = sum(sensN_bR, 2);
    sens_NT(1:nt, 3, iexp)  = sum(sensN_alpha_b, 2);
    sens_NT(1:nt, 4, iexp)  = sum(sensN_d_maxS, 2);
    sens_NT(1:nt, 5, iexp)  = sum(sensN_alpha_d, 2);
    sens_NT(1:nt, 6, iexp)  = sum(sensN_beta_d, 2);
    sens_NT(1:nt, 7, iexp)  = sum(sensN_EC_50d, 2);
    sens_NT(1:nt, 8, iexp)  = sum(sensN_H_d, 2);
    sens_NT(1:nt, 9, iexp)  = sum(sensN_xiSR, 2);
    sens_NT(1:nt, 10, iexp) = sum(sensN_kxi, 2);
    sens_NT(1:nt, 11, iexp) = sum(sensN_NT0, 2);
    sens_NT(1:nt, 12, iexp) = sum(sensN_lambda_T0, 2);
    
end

          
end