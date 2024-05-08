%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System of ODEs for state and approximate variance (constant control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dsdt = Odes_Var(tt, ss, tmod, xx, mug, muk, munT)
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

munT = interp1(tmod, munT, tt);
xx   = interp1(tmod, xx, tt);

% --------------------------------------------------------------------- %
% ODEs for variance of the approximate process:
varxT = ss;                                                        % Variance of xT,

% ODE for variance of xT:
dsdt  = (mug + muk).'*xx.' + 2*munT*varxT;
    


end