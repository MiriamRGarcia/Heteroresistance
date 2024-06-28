%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Odes_aveBD_Sens: System of ODEs for state sensitivities of the
%                  average BD model with 
%                  constant antimicrobial concentration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dsdt = Odes_aveBD_Sens(t, s, AA, AA_bS, AA_bR, AA_alpha_b, AA_dmaxS,...
                AA_alpha_d, AA_beta_d, AA_EC50d, AA_Hd, AA_xiSR, AA_kxi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT :
% t        = Time variable,
% s        = State variables, s = (N, N_theta=dN/dtheta),
% AA       = Coefficient matrix,
% AA_theta = Derivative of coefficient matrix AA with respect to theta,
%
% OUTPUT:
% dsdt = Time derivatives of state sensitivities;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem sizes:
nr = size(AA, 2);

% ----------------------------------------------------------------------- %
% Odes for the average process:
N    = s(1:nr);                                                    
dsdt = AA*N;                                                             

% ----------------------------------------------------------------------- %
% Odes for sensitivities of the average proces:

% Derivative with respect to bS:
sbS      = s((nr + 1):2*nr);
dsdt     = [dsdt;AA_bS*N + AA*sbS];

% Derivative with respect to mugR:
sbR      = s((2*nr + 1):3*nr);
dsdt     = [dsdt;AA_bR*N + AA*sbR];

% Derivative with respect to alph_g:
salpha_b = s((3*nr + 1:4*nr));
dsdt     = [dsdt;AA_alpha_b*N + AA*salpha_b];

% Derivative with respect to mukmaxS:
sdmaxS   = s((4*nr + 1):5*nr);
dsdt     = [dsdt;AA_dmaxS*N + AA*sdmaxS];

% Derivative with respect to alpha_d:
salpha_d  = s((5*nr + 1):6*nr);
dsdt     = [dsdt; AA_alpha_d*N + AA*salpha_d];

% Derivative with respect to beta_d:
sbeta_d  = s((6*nr + 1):7*nr);
dsdt     = [dsdt; AA_beta_d*N + AA*sbeta_d];

% Derivative with respect to EC50d:
sEC50d   = s((7*nr + 1):8*nr);
dsdt     = [dsdt; AA_EC50d*N + AA*sEC50d];

% Derivative with respect to Hk:
sHd      = s((8*nr + 1):9*nr);
dsdt     = [dsdt; AA_Hd*N + AA*sHd];

% Derivative with respect to xiSR:
sxi_SR   = s((9*nr + 1):10*nr);
dsdt     = [dsdt;AA_xiSR*N + AA*sxi_SR];

% Derivative with respect to kxi:
sk_xi    = s((10*nr + 1):11*nr);
dsdt     = [dsdt;AA_kxi*N + AA*sk_xi];

% Derivative with respect N_T0:
sN_T0    = s(11*nr + 1:12*nr);
dsdt     = [dsdt;AA*sN_T0];

% Derivative with respect lambda_T0:
slambda_T0 = s(12*nr + 1:13*nr);
dsdt       = [dsdt;AA*slambda_T0];


end