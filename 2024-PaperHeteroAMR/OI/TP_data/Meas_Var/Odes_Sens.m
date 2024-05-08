%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System of ODEs for state and variance sensitivities (constant control)
% for the approximate process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dsdt = Odes_Sens(tt, ss, AA, AA_mugS0, AA_mugR, ...
                AA_alph_g, AA_mukmaxS, AA_r50k, AA_alph_k, AA_EC50k,...
                AA_Hk, AA_xiSR, AA_kxi, mug, muk, der_mug, der_muk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT :
% tt       = Time variable,
% ss       = State variables, ss = (xx, xx_theta, var, var_theta),
% AA       = Coefficient matrix,
% AA_theta = Derivative of coefficient matrix AA with respect to theta,
% mug      = Array of growth rates,
% muk      = Array of kill rates,
% der_mug  = Array nr x np with derivatives of mug with respect theta,
% der_muk  = Array nr x np with derivatives of muk with respect theta,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of AMR degrees:
nr   = size(AA, 2);

% ----------------------------------------------------------------------- %
% Odes for average process
% ----------------------------------------------------------------------- %
xx   = ss(1:nr);                                                           % Array with averages for the AMR degrees,
dsdt = AA*xx;                                                              % ODEs for averages,

% ----------------------------------------------------------------------- %
% Odes for sensitivities of the average proces
% ----------------------------------------------------------------------- %

% Derivative with respect to mugS:
smugS    = ss((nr + 1):2*nr);
dsdt     = [dsdt;AA_mugS0*xx + AA*smugS];

% Derivative with respect to mugR:
smugR    = ss((2*nr + 1):3*nr);
dsdt     = [dsdt;AA_mugR*xx + AA*smugR];

% Derivative with respect to alph_g:
salph_g  = ss((3*nr + 1:4*nr));
dsdt     = [dsdt;AA_alph_g*xx + AA*salph_g];

% Derivative with respect to mukmaxS:
smukmaxS = ss((4*nr + 1):5*nr);
dsdt     = [dsdt;AA_mukmaxS*xx + AA*smukmaxS];

% Derivative with respect to r50k:
sr50k    = ss((5*nr + 1):6*nr);
dsdt     = [dsdt; AA_r50k*xx + AA*sr50k];

% Derivative with respect to alph_k:
salph_k  = ss((6*nr + 1):7*nr);
dsdt     = [dsdt; AA_alph_k*xx + AA*salph_k];

% Derivative with respect to EC50k:
sEC50k   = ss((7*nr + 1):8*nr);
dsdt     = [dsdt; AA_EC50k*xx + AA*sEC50k];

% Derivative with respect to Hk:
sHk      = ss((8*nr + 1):9*nr);
dsdt     = [dsdt; AA_Hk*xx + AA*sHk];

% Derivative with respect to xiSR:
sxiSR    = ss((9*nr + 1):10*nr);
dsdt     = [dsdt;AA_xiSR*xx + AA*sxiSR];

% Derivative with respect to kxi:
skxi     = ss((10*nr + 1):11*nr);
dsdt     = [dsdt;AA_kxi*xx + AA*skxi];

% Derivative with respect X0:
sX0      = ss(11*nr + 1:12*nr);
dsdt     = [dsdt;AA*sX0];

% Derivative with respect lamb_IC:
slamb_IC = ss(12*nr + 1:13*nr);
dsdt     = [dsdt;AA*slamb_IC];

% ----------------------------------------------------------------------- %
% ODEs for variance of the approximate process
% ----------------------------------------------------------------------- %
xT    = sum(xx);                                                           % Total number of cells,
munT  = (mug - muk).'*xx/xT;                                               % Total net growth rate,
varxT = ss(13*nr + 1);                                                     % Variance of xT,

% ODE for variance of xT:
dsdt  = [dsdt;(mug + muk).'*xx + 2*munT*varxT];

% ----------------------------------------------------------------------- %
% ODEs for sensitivities of variance for the approximate process
% ----------------------------------------------------------------------- %
svarxT    = ss(13*nr + 2:13*(nr + 1));                                     % Array of variance sensitivities to parameters,

sst       = [smugS smugR salph_g smukmaxS sr50k salph_k sEC50k ...
             sHk sxiSR skxi sX0 slamb_IC];                                 % Matrix with state sensitivities
sxT       = sum(sst.', 2);

sens_munT = ((der_mug - der_muk).'*xx + sst.'*(mug - muk))/xT - (mug - muk).'*xx*sxT/(xT^2);

dsvardt   = 2*varxT*sens_munT + 2*munT*svarxT + (der_mug + der_muk).'*xx +  sst.'*(mug + muk);


dsdt = [dsdt;
        dsvardt];

end