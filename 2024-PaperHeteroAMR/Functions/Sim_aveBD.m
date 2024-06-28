%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sim_aveBD: Simulate average BD heteroresistance model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N, N_T] = Sim_aveBD(r, R, tsim, Cexp, pars, ODEoptions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% r           = Discretisation of the AMR level (values between entire 
%               sensitivity 0 and entire resistance 1) (size: m_r x 1),
% R           = Auxiliary matrix with jumps in AMR level (size: m_r x m_r);
% tsim        = Array of simulation times for solving ODEs (size: m_t x 1);
% Cexp        = Array of antimicrobial concentrations (assumed  
%               constant at each experiment) (size: m_e x 1),
% pars       = Values of parameters used for simulation (size: m_p x 1);
% ODEoptions = Options for ODE solver (ode15s);
%
% OUTPUT:
% N_T        = Total average counts at the different experiments 
%             (size: m_t x m_e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem sizes:
m_t = numel(tsim);
m_e = numel(Cexp);
m_r = numel(r);

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
    
% Calculate initial condition:  
f_0 = exp(-lambda_T0*r);
f_0 = f_0/sum(f_0);
N_0 = N_T0*f_0;                                                 
   
% Calculate auxiliary coefficient matrix:                                                                                                       
Xi     = xi_SR*exp(k_xi*(1 - R));
Xi     = Xi - diag(diag(Xi));                                  
AA_aux = Xi.' - diag(sum(Xi, 2));

% Calculate birth and maximal kill rates:
b     = b_S*b_R./(b_R + r.^alpha_b*(b_S - b_R));                
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);                                  

% Calculate cell counts in each experiment:
N_T = zeros(m_t, m_e);
N   = zeros(m_t, m_r, m_e);

for iexp = 1:m_e
    
    % Calculate kill rate:
    C  = Cexp(iexp);  
    HC = C^H_d/(C^H_d + EC_50d^H_d);
    d  = d_max*HC;                                             
    
    % Coefficient matrix of linear ODEs:
    AA = AA_aux + diag(b - d);                                    
    
    % Obtain cell counts at the current experiment:
    [~, xout] = ode15s(@(t,s) Odes_aveBD(t, s, AA), tsim, N_0, ODEoptions);
    
    % Almacenate cell counts at the current experiment:
    N(1:m_t, 1:m_r, iexp) = xout;
    
    % Almacenate total cell counts at the current experiment:
    N_T(1:m_t, iexp)      = sum(xout, 2);                   
end 
        