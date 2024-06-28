%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sim_aveBD_2subpop: Simulate average BD heteroresistance model with
%                    only two subpopulations (m_r = 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N, N_T] = Sim_aveBD_2subpop(r, R, tsim, Cexp, pars, ODEoptions)
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
m_r       = numel(r);
m_t       = numel(tsim);
m_e       = numel(Cexp);

% Obtain parameter values:
b_S       = pars(1);                                                      
b_R       = pars(2);                                                     

d_maxS    = pars(3);                                                    
EC_50d    = pars(4);                                                        
H_d       = pars(5);                                                           

xi_SR     = pars(6);     

N_T0      = pars(7);                                                       
lambda_T0 = pars(8);   
    
% Calculate initial condition:  
f_0 = exp(-lambda_T0*r);
f_0 = f_0/sum(f_0);
N_0 = N_T0*f_0;                                                 
   
% Calculate auxiliary coefficient matrix:                                                                                                       
Xi     = [0 xi_SR;xi_SR 0];                                        
Xi     = Xi - diag(diag(Xi));
AA_aux = Xi.' - diag(sum(Xi, 2));                                  

% Calculate birth rate:
b     = [b_S;b_R];                

% Calculate cell counts in each experiment:
N_T = zeros(m_t, m_e);
N   = zeros(m_t, m_r, m_e);

for iexp = 1:m_e
    
    % Calculate kill rate:
    C  = Cexp(iexp);  
    HC = C^H_d/(C^H_d + EC_50d^H_d);
    d  = [d_maxS*HC;0];                                            
    
    % Coefficient matrix of linear ODEs:
    AA = AA_aux + diag(b - d);                                    
    
    % Obtain cell counts at the current experiment:
    [~, xout] = ode15s(@(t,s) Odes_aveBD(t, s, AA), tsim, N_0, ODEoptions);
    
    % Almacenate cell counts at the current experiment:
    N(1:m_t, 1:m_r, iexp) = xout;
    
    % Almacenate total cell counts at the current experiment:
    N_T(1:m_t, iexp) = sum(xout, 2);                   
end 
        