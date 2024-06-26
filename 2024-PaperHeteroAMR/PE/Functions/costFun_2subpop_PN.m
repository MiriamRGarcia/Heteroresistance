%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% costFun_2subpop_MNHe: Cost function for MLE in PN case (2 subpopulations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = costFun_2subpop_PN(decVars, r, R, tmod, texp_ind, Cexp,...
                      N_Tave_data, pars_nom, Var_data, Weights, ODEoptions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% decVars    = Decision variables (parameters to calibrate);
% r          = Discretisation of the AMR level (values between entire 
%            sensitivity S = 0 and entire resistance R = 1) (nr x 1),
% tmod       = Simulation times for solve ODEs (nt x 1);
% texp_ind   = Indexes of sampling times (ntexp x 1);
% Cexp       = Antimicrobial concentrations (assumed constant 
%            at each experiment),
% N_Tdata    = Total count data at sampling times for MLE,
% pars_nom   = Nomimal values of parameters;
% Var_data   = Variance in total count data at sampling times (for PN case);
% Weights    = Weights to remove NaN data (for PN case);
% ODEoptions = Options for ODE solver;
%
% OUTPUT:
% J = Value of the cost function;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% Obtain parameter values:
if size(pars_nom) == size(decVars)
    decVars = pars_nom.*decVars;
else
    decVars = pars_nom.'.*decVars;
end

b_S       = decVars(1);
b_R       = decVars(2);
d_maxS    = decVars(3);
EC_50d    = decVars(4);
H_d       = decVars(5);
xi_SR     = decVars(6);
N_T0      = decVars(7);
lambda_T0 = decVars(8);

% ----------------------------------------------------------------------- %
% Previous calculations:
f0  = exp(-lambda_T0*r);
f0  = f0/sum(f0);
N_0 = N_T0*f0;                                                             

Xi     = [xi_SR 0;0 xi_SR];

AA_aux = Xi.' - diag(sum(Xi, 2));

b      = [b_S;b_R];
d_max  = [d_maxS;0];

% ----------------------------------------------------------------------- %
% Solve ODEs:
Nexp   = numel(Cexp);
ntexp  = numel(texp_ind);

N_Tmod = zeros(ntexp, Nexp);

for iexp = 1:Nexp
    
    C  = Cexp(iexp);
    d  = d_max*C^H_d/(C^H_d + EC_50d^H_d);

    AA = AA_aux + diag(b - d);
    
    [~, xout] = ode15s(@(s,y) Odes_cte(s, y, AA), tmod, N_0, ODEoptions);
    
    N_T                   = sum(xout, 2);
    N_Tmod(1:ntexp, iexp) = N_T(texp_ind);

end

% ----------------------------------------------------------------------- %
% Calculate cost function:
J  = sum(sum(Weights.*((N_Tave_data - N_Tmod).^2)./Var_data));

end
