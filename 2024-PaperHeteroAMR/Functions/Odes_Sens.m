%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Odes_Sens: System of ODEs for the state sensitivities of the
%            average BD heteroresistance model with 
%            constant antimicrobial concentration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dsdt = Odes_Sens(t, s, AA_Sens)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT :
% t        = Time variable,
% s        = State variables, s = (N, N_theta=dN/dtheta),
% AA       = Array of size m_r x m_r x m_p + 1 with coefficient matrix
%            and matrix derivatives with respect to parameters;
% OUTPUT:
% dsdt = Time derivatives of state sensitivities;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem sizes:
m_r     = size(AA_Sens, 2);
m_p_aux = size(AA_Sens, 3) - 1;
m_p     = m_p_aux + 2;

% ----------------------------------------------------------------------- %
% Odes for the average of cell counts:
N    = s(1:m_r);     
AA   = AA_Sens(1:m_r, 1:m_r, 1);

dsdt = AA*N;     

% ----------------------------------------------------------------------- %
% Odes for sensitivities of the average cell counts:

% Loop in the model parameters:
for ip = 1:m_p_aux   
    sens_aux    = s((ip*m_r + 1):(ip + 1)*m_r);
    AA_sens_aux = reshape(AA_Sens(1:m_r, 1:m_r, 1 + ip), m_r, m_r);
    dsdt        = [dsdt;AA_sens_aux*N + AA*sens_aux];
end

% Ampliate with state sensitivities to N_T0:
sN_T0    = s((m_p - 1)*m_r + 1:m_p*m_r);
dsdt     = [dsdt;AA*sN_T0];

% Ampliate with state sensitivities to lambda_T0:
slambda_T0 = s(m_p*m_r + 1:(m_p + 1)*m_r);
dsdt       = [dsdt;AA*slambda_T0];


end