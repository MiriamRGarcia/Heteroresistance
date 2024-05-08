%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System of ODEs for state and approximate variance (constant control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dsdt = Odes_MeanVar(tt, ss, AA, mug, muk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT :
% tt       = Time variable,
% ss       = State variables, ss = (xx, varxxT),
% AA       = Coefficient matrix,
% mug      = Array of growth rates,
% muk      = Array of kill rates,
% der_mug  = Array nr x np with derivatives of mug with respect theta,
% der_muk  = Array nr x np with derivatives of muk with respect theta,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% Number of AMR degrees:
nr = size(AA, 2);

% ----------------------------------------------------------------------- %
% Odes for average process:
xx   = ss(1:nr);                                                           % Array with averages for the AMR degrees,
dsdt = AA*xx;                                                              % ODEs for averages,

% ----------------------------------------------------------------------- %
% ODEs for variance of the approximate process:
xT    = sum(xx);                                                           % Total number of cells,
munT  = (mug - muk).'*xx/xT;                                               % Total net growth rate,
varxT = ss(nr + 1);                                                        % Variance of xT,

% ODE for variance of xT:
dsdt  = [dsdt;
        (mug + muk).'*xx + 2*munT*varxT];
    


end