%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System of ODEs for state and variance sensitivities (constant control)
% for the approximate process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dsdt = Odes_MV(tt, ss, AA, mug, muk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT :
% tt       = Time variable,
% ss       = State variables, ss = (xx, xx_theta, var, var_theta),
% AA       = Coefficient matrix,
% mug      = Array of growth rates,
% muk      = Array of kill rates,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of AMR degrees:
nr   = size(AA, 2);

% ----------------------------------------------------------------------- %
% Odes for average process
% ----------------------------------------------------------------------- %
xx   = ss(1:nr);                                                           % Array with averages for the AMR degrees,
dsdt = AA*xx;                                                              % ODEs for averages,


% ----------------------------------------------------------------------- %
% ODEs for variance of the approximate process
% ----------------------------------------------------------------------- %
xT    = sum(xx);                                                           % Total number of cells,
munT  = (mug - muk).'*xx/xT;                                               % Total net growth rate,
varxT = ss(nr + 1);                                                        % Variance of xT,

% ODE for variance of xT:
dsdt  = [dsdt;
        (mug + muk).'*xx + 2*munT*varxT];


end