%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SensME : Function calculating the state sensitivities of the average
%          BD heteroresistance model for  
%          multiple experiments at (constant) antimicrobial concentrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NT, sens_NT] = SensME(tsim, r, Cexp, pars, N_0, ODEoptions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% tsim       = Simulation times for solving ODEs (m_t x 1);
% r          = Discretisation of the AMR level (values between entire 
%              sensitivity = 0 and entire resistance = 1) (size: m_r x 1),
% Cexp       = Array of constant antimicrobial concentrations (size: m_e x 1),
% pars       = Parameter values;
% N_0        = Initial cell counts (size: m_r x 1);
% ODEoptions = Options for ODEs solver (ode15s);
%
% OUTPUT:
% N_T     = Total average counts at simulation times (m_t x m_e);
% sens_NT = State sensitivities of total average counts (size: m_t x m_p x m_e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% Problem sizes:                                                     
m_t = numel(tsim);                                                    
m_r = numel(r);
m_e = numel(Cexp);
m_p = numel(pars);                                                   

% ----------------------------------------------------------------------- %
% Obtain parameter values and previous calculations:

% Auxiliary matrix to calculate mutation rates:
RR  = repmat(r, 1, m_r) - repmat(r.', m_r, 1);                         
RR  = RR - triu(RR) + tril(RR).';

% Initial condition for average cell counts and state sensitivities:
s0 = [N_0;
      zeros(m_r*(m_p - 2),1)];                                             % Initial conditions are independent of parameters
                                                                           % unless for X0 and lamb_IC
% Ampliate with the initial condition for sensitivity to IC parameters:
N_T0       = pars(m_p - 1);
lambda_T0  = pars(m_p);

f0_aux     = exp(-lambda_T0*r);
s0_N0      = f0_aux/sum(f0_aux);
RR_aux     = repmat(r.', m_r, 1) - repmat(r, 1, m_r);    
s0_lamb_IC = N_T0/(sum(f0_aux)^2)*f0_aux.*(RR_aux*f0_aux);

s0 = [s0;
      s0_N0;
      s0_lamb_IC];

% Initialise total cell counts and state sensitivities:
NT         = zeros(m_t, m_e);
sens_NT    = zeros(m_t, m_p, m_e);

% Initialise derivatives with respect common parameters:
d_d_maxS   = zeros(m_r,1);
d_EC_50d   = zeros(m_r,1);
d_H_d      = zeros(m_r,1);

if m_r > 2
    
    % Parameter values:
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
    
    % Mutation rates matrix:
    Xi = xi_SR*exp(k_xi*(1 - RR));
    Xi = Xi - diag(diag(Xi));

    % Auxiliary coefficient matrix for average cell counts:
    AA_aux = Xi.' - diag(sum(Xi, 2));

    % Derivative of the dynamics with respect to bS:
    d_bS(1)         = 1;
    d_bS(m_r)       = 0;
    d_bS(2:(m_r-1)) = bR^2*(1 - r(2:(m_r-1)).^alpha_b)./(r(2:(m_r-1)).^alpha_b*(bR - bS) - bR).^2;

    AA_bS = diag(d_bS);

    % Derivative with respect to bR:
    d_bR(1)         = 0;
    d_bR(m_r)       = 1;
    d_bR(2:(m_r-1)) = bS^2*r(2:(m_r-1)).^alpha_b./(bR + r(2:(m_r-1)).^alpha_b*(bS - bR)).^2;

    AA_bR = diag(d_bR);

    % Derivative with respect to alpha_b:
    d_alpha_b(1)         = 0;
    d_alpha_b(m_r)       = 0;
    d_alpha_b(2:(m_r-1)) = bR*bS*(bR - bS)*log(r(2:(m_r-1))).*r(2:(m_r-1)).^alpha_b./(bR + r(2:(m_r-1)).^alpha_b*(bS - bR)).^2;

    AA_alpha_b  = diag(d_alpha_b);

    % Derivative with respect to xi_SR:
    d_xi_SR  = exp(k_xi*(1 - RR));
    d_xi_SR  = d_xi_SR - diag(diag(d_xi_SR));
    
    AA_xi_SR = d_xi_SR - diag(sum(Xi/xi_SR, 2));

    % Derivative with respect to k_xi:
    d_k_xi  = (1 - RR).*Xi;
    d_k_xi  = d_k_xi - diag(diag(d_k_xi));
    
    AA_k_xi = d_k_xi - diag(sum(d_k_xi, 2));

    % Calculate birth rate:
    b     = bS*bR./(bR + (bS - bR)*r.^alpha_b);
 
    % Calculate maximal kill rate:
    d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);    
    
    % Calculate derivatives of coefficient matrix with respect to
    % parameters varying with antimicrobial concentration:
    d_alpha_d  = zeros(m_r,1);
    d_beta_d   = zeros(m_r,1);

    for iexp = 1:m_e
    
        %-------------------------------------------------%
        % Calculate coefficient matrix of the state system:
    
        % Drug concentration:
        CC = Cexp(iexp);
    
        % Calculate kill rate at current time:
        HC = CC^H_d/(CC^H_d + EC_50d^H_d);
        d  = d_max*HC;

        AA = AA_aux + diag(b - d);
    
        %-------------------------------------------------%
        % Calculate derivatives of the coefficient matrix:
    
        % Derivative of A with respect to d_maxS:
        d_d_maxS(1)   = 1;
        d_d_maxS(m_r) = 0;
        d_d_maxS(2:(m_r - 1)) = beta_d^alpha_d*(1 - r(2:(m_r - 1)).^alpha_d)./(beta_d^alpha_d + r(2:(m_r - 1)).^alpha_d);

        d_d_maxS  = - HC*d_d_maxS;

        AA_d_maxS = diag(d_d_maxS);
    
        % Derivative with respect to alpha_d:
        d_alpha_d(1)   = 0;
        d_alpha_d(m_r) = 0;
        
        d_alpha_d(2:(m_r -1)) = beta_d^alpha_d*d_maxS*r(2:(m_r -1)).^alpha_d.*((1 + beta_d^alpha_d)*log(r(2:(m_r - 1))) + (r(2:(m_r -1)).^alpha_d - 1)*log(beta_d))./(r(2:(m_r -1)).^alpha_d + beta_d^alpha_d).^2;

        d_alpha_d  = HC*d_alpha_d;

        AA_alpha_d = diag(d_alpha_d);
    
        % Derivative with respect to beta_d:
        d_beta_d(1)          = 0;
        d_beta_d(m_r)        = 0;
        d_beta_d(2:(m_r -1)) =  alpha_d*beta_d^(alpha_d - 1)*d_maxS*(r(2:(m_r -1)).^alpha_d - 1).*r(2:(m_r -1)).^alpha_d./(beta_d^alpha_d + r(2:(m_r -1)).^alpha_d).^2;

        d_beta_d  = HC*d_beta_d;
    
        AA_beta_d = diag(d_beta_d);

        % Derivative with respect to EC50_d:
        d_EC_50d(m_r)         = 0;
        d_EC_50d(1:(m_r - 1)) = - EC_50d^(H_d - 1)*H_d*beta_d^alpha_d*d_maxS*(r(1:(m_r -1)).^alpha_d - 1)./(r(1:(m_r -1)).^alpha_d + beta_d^alpha_d);
        d_EC_50d              = HC^2*d_EC_50d;

        AA_EC_50d = diag(d_EC_50d);

        % Derivative with respect to H_d:
        if CC > 1e-10
            d_H_d(1:(m_r -1)) = EC_50d^H_d*beta_d^alpha_d*d_maxS*(log(CC) - log(EC_50d))*(r(1:(m_r -1)).^alpha_d - 1)./(r(1:(m_r -1)).^alpha_d + beta_d^alpha_d);
            d_H_d(m_r)        = 0;
            d_H_d             = HC^2*d_H_d;
        else
            d_H_d(1:m_r) = 0;
        end
    
        AA_H_d               = diag(d_H_d);
    
        %-------------------------------------------------%
        % Calculate output sensitivities:
        
        % Tensor with derivatives of coefficient matrix:
        AA_Sens(1:m_r,1:m_r, 1)  = AA;
        AA_Sens(1:m_r,1:m_r, 2)  = AA_bS;
        AA_Sens(1:m_r,1:m_r, 3)  = AA_bR;
        AA_Sens(1:m_r,1:m_r, 4)  = AA_alpha_b;
        AA_Sens(1:m_r,1:m_r, 5)  = AA_d_maxS;
        AA_Sens(1:m_r,1:m_r, 6)  = AA_alpha_d;
        AA_Sens(1:m_r,1:m_r, 7)  = AA_beta_d;
        AA_Sens(1:m_r,1:m_r, 8)  = AA_EC_50d;
        AA_Sens(1:m_r,1:m_r, 9)  = AA_H_d;
        AA_Sens(1:m_r,1:m_r, 10) = AA_xi_SR;
        AA_Sens(1:m_r,1:m_r, 11) = AA_k_xi;
            
        % Call to ODEs giving state sensitivities for the current
        % experiment:
        [~, St_sens] = ode15s(@(t,s) Odes_Sens(t, s, AA_Sens), tsim, s0, ODEoptions);
   
        % Total population:
        N               = St_sens(1:m_t, 1:m_r);
        NT(1:m_t, iexp) = sum(N, 2);
      
        % State sensitivities:
        sensN_bS        = St_sens(1:m_t,m_r + 1:2*m_r);
        sensN_bR        = St_sens(1:m_t,2*m_r + 1:3*m_r);
        sensN_alpha_b   = St_sens(1:m_t,3*m_r + 1:4*m_r);
        sensN_d_maxS    = St_sens(1:m_t,4*m_r + 1:5*m_r);
        sensN_alpha_d   = St_sens(1:m_t,5*m_r + 1:6*m_r);
        sensN_beta_d    = St_sens(1:m_t,6*m_r + 1:7*m_r);
        sensN_EC_50d    = St_sens(1:m_t,7*m_r + 1:8*m_r);
        sensN_H_d       = St_sens(1:m_t,8*m_r + 1:9*m_r);
        sensN_xiSR      = St_sens(1:m_t,9*m_r + 1:10*m_r);
        sensN_kxi       = St_sens(1:m_t,10*m_r + 1:11*m_r);
        sensN_NT0       = St_sens(1:m_t,11*m_r + 1:12*m_r);
        sensN_lambda_T0 = St_sens(1:m_t,12*m_r + 1:13*m_r);
    
        % Sensitivity of the total population:
        sens_NT(1:m_t, 1, iexp)  = sum(sensN_bS, 2);
        sens_NT(1:m_t, 2, iexp)  = sum(sensN_bR, 2);
        sens_NT(1:m_t, 3, iexp)  = sum(sensN_alpha_b, 2);
        sens_NT(1:m_t, 4, iexp)  = sum(sensN_d_maxS, 2);
        sens_NT(1:m_t, 5, iexp)  = sum(sensN_alpha_d, 2);
        sens_NT(1:m_t, 6, iexp)  = sum(sensN_beta_d, 2);
        sens_NT(1:m_t, 7, iexp)  = sum(sensN_EC_50d, 2);
        sens_NT(1:m_t, 8, iexp)  = sum(sensN_H_d, 2);
        sens_NT(1:m_t, 9, iexp)  = sum(sensN_xiSR, 2);
        sens_NT(1:m_t, 10, iexp) = sum(sensN_kxi, 2);
        sens_NT(1:m_t, 11, iexp) = sum(sensN_NT0, 2);
        sens_NT(1:m_t, 12, iexp) = sum(sensN_lambda_T0, 2);
    
    end    
else
    bS        = pars(1);
    bR        = pars(2);

    d_maxS    = pars(3);
    EC_50d    = pars(4);
    H_d       = pars(5);

    xi_SR     = pars(6);
    
    % Mutation rates matrix:
    Xi = [0 xi_SR;xi_SR 0];
    Xi = Xi - diag(diag(Xi));

    % Auxiliary coefficient matrix for average cell counts:
    AA_aux = Xi.' - diag(sum(Xi, 2));

    % Derivative of the dynamics with respect to bS:
    d_bS(1)         = 1;
    d_bS(m_r)       = 0;

    AA_bS = diag(d_bS);

    % Derivative with respect to bR:
    d_bR(1)         = 0;
    d_bR(m_r)       = 1;

    AA_bR = diag(d_bR);
    
    % Derivative with respect to xi_SR:
    AA_xi_SR = [-1 1;1 -1];

    % Calculate birth rate:
    b     = [bS;bR];
 
    % Calculate maximal kill rate:
    d_max = [d_maxS;0]; 
    
    for iexp = 1:m_e
    
        %-------------------------------------------------%
        % Calculate coefficient matrix of the state system:
    
        % Drug concentration:
        CC = Cexp(iexp);
    
        % Calculate kill rate at current time:
        HC = CC^H_d/(CC^H_d + EC_50d^H_d);
        d  = d_max*HC;

        AA = AA_aux + diag(b - d);
    
        %-------------------------------------------------%
        % Calculate derivatives of the coefficient matrix:
        
        % Derivative with respect to d_maxS:
        d_d_maxS(1)   = - HC;
        d_d_maxS(m_r) = 0;
    
        AA_d_maxS = diag(d_d_maxS);
    
        % Derivative with respect to EC50_d:
        d_EC_50d(1)   = EC_50d^(H_d - 1)*H_d*d_maxS/(CC^H_d + EC_50d^H_d)*HC;
        d_EC_50d(m_r) = 0;
        
        AA_EC_50d = diag(d_EC_50d);

        % Derivative with respect to H_d:
        d_H_d(m_r) = 0;
        if CC > 1e-10
            d_H_d(1) = - EC_50d^H_d*d_maxS*(log(CC) - log(EC_50d))/(CC^H_d + EC_50d^H_d)*HC;
        else
            d_H_d(1) = 0;
        end
    
        AA_H_d = diag(d_H_d);
    
        %-------------------------------------------------%
        % Calculate output sensitivities:
        
        % Tensor with derivatives of coefficient matrix:
        AA_Sens(1:m_r,1:m_r, 1) = AA;
        AA_Sens(1:m_r,1:m_r, 2) = AA_bS;
        AA_Sens(1:m_r,1:m_r, 3) = AA_bR;
        AA_Sens(1:m_r,1:m_r, 4) = AA_d_maxS;
        AA_Sens(1:m_r,1:m_r, 5) = AA_EC_50d;
        AA_Sens(1:m_r,1:m_r, 6) = AA_H_d;
        AA_Sens(1:m_r,1:m_r, 7) = AA_xi_SR;
        
        
        % Call to ODEs giving state sensitivities for the current
        % experiment:
        [~, St_sens] = ode15s(@(t,s) Odes_Sens(t, s, AA_Sens), tsim, s0, ODEoptions);
   
        % Total population:
        N               = St_sens(1:m_t, 1:m_r);
        NT(1:m_t, iexp) = sum(N, 2);
      
        % State sensitivities:
        sensN_bS        = St_sens(1:m_t,m_r + 1:2*m_r);
        sensN_bR        = St_sens(1:m_t,2*m_r + 1:3*m_r);
        sensN_d_maxS    = St_sens(1:m_t,3*m_r + 1:4*m_r);
        sensN_EC_50d    = St_sens(1:m_t,4*m_r + 1:5*m_r);
        sensN_H_d       = St_sens(1:m_t,5*m_r + 1:6*m_r);
        sensN_xiSR      = St_sens(1:m_t,6*m_r + 1:7*m_r);
        sensN_NT0       = St_sens(1:m_t,7*m_r + 1:8*m_r);
        sensN_lambda_T0 = St_sens(1:m_t,8*m_r + 1:9*m_r);
    
        % Sensitivity of the total population:
        sens_NT(1:m_t, 1, iexp) = sum(sensN_bS, 2);
        sens_NT(1:m_t, 2, iexp) = sum(sensN_bR, 2);
        sens_NT(1:m_t, 3, iexp) = sum(sensN_d_maxS, 2);
        sens_NT(1:m_t, 4, iexp) = sum(sensN_EC_50d, 2);
        sens_NT(1:m_t, 5, iexp) = sum(sensN_H_d, 2);
        sens_NT(1:m_t, 6, iexp) = sum(sensN_xiSR, 2);
        sens_NT(1:m_t, 7, iexp) = sum(sensN_NT0, 2);
        sens_NT(1:m_t, 8, iexp) = sum(sensN_lambda_T0, 2);
    
    end 
    
end

          
end