%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OdesSens: System of ODEs for state sensitivities (constant control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dsdt = OdesSens(t, ss, AA, AA_bS, AA_bR, AA_alpha_b, AA_dmaxS,...
                AA_alpha_d, AA_beta_d, AA_EC50d, AA_Hd, AA_xiSR, AA_kxi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT :
% t        = Time variable,
% ss       = State variables, ss = (N, N_theta),
% AA       = Coefficient matrix,
% AA_theta = Derivative of coefficient matrix AA with respect to theta,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of AMR levels:
nr   = size(AA, 2);

% ----------------------------------------------------------------------- %
% Odes for the average process:
N    = ss(1:nr);                                                    
dsdt = AA*N;                                                             

% ----------------------------------------------------------------------- %
% Odes for sensitivities of the average proces:

% Derivative with respect to bS:
sb_S     = ss((nr + 1):2*nr);
dsdt     = [dsdt;AA_bS*N + AA*sb_S];

% Derivative with respect to mugR:
sb_R     = ss((2*nr + 1):3*nr);
dsdt     = [dsdt;AA_bR*N + AA*sb_R];

% Derivative with respect to alph_g:
salpha_b = ss((3*nr + 1:4*nr));
dsdt     = [dsdt;AA_alpha_b*N + AA*salpha_b];

% Derivative with respect to mukmaxS:
sdmax_S  = ss((4*nr + 1):5*nr);
dsdt     = [dsdt;AA_dmaxS*N + AA*sdmax_S];

% Derivative with respect to alpha_d:
salpha_d  = ss((5*nr + 1):6*nr);
dsdt     = [dsdt; AA_alpha_d*N + AA*salpha_d];

% Derivative with respect to beta_d:
sbeta_d  = ss((6*nr + 1):7*nr);
dsdt     = [dsdt; AA_beta_d*N + AA*sbeta_d];

% Derivative with respect to EC50d:
sEC50d   = ss((7*nr + 1):8*nr);
dsdt     = [dsdt; AA_EC50d*N + AA*sEC50d];

% Derivative with respect to Hk:
sHd      = ss((8*nr + 1):9*nr);
dsdt     = [dsdt; AA_Hd*N + AA*sHd];

% Derivative with respect to xiSR:
sxi_SR   = ss((9*nr + 1):10*nr);
dsdt     = [dsdt;AA_xiSR*N + AA*sxi_SR];

% Derivative with respect to kxi:
sk_xi    = ss((10*nr + 1):11*nr);
dsdt     = [dsdt;AA_kxi*N + AA*sk_xi];

% Derivative with respect N_T0:
sN_T0    = ss(11*nr + 1:12*nr);
dsdt     = [dsdt;AA*sN_T0];

% Derivative with respect lambda_T0:
slambda_T0 = ss(12*nr + 1:13*nr);
dsdt       = [dsdt;AA*slambda_T0];


end