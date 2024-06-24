%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% costFun_MNHe: Cost function for MLE in MNHe case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = costFun_MNHe(decVars, r, R, tmod, texp_ind, Cexp, N_Tave_data, pars_nom, Var_data, Weights, ODEoptions)
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

% Obtain parameter values:
if size(pars_nom) == size(decVars)
    decVars = pars_nom.*decVars;
else
    decVars = pars_nom.'.*decVars;
end

b_S       = decVars(1);
b_R       = decVars(2);
alpha_b   = decVars(3);
d_maxS    = decVars(4);
alpha_d   = decVars(5);
beta_d    = decVars(6);
EC_50d    = decVars(7);
H_d       = decVars(8);
xi_SR     = decVars(9);
k_xi      = decVars(10);
N_T0      = decVars(11);
lambda_T0 = decVars(12);
var_b     = decVars(13);

% ----------------------------------------------------------------------- %
% Previous calculations:
f0  = exp(-lambda_T0*r);
f0  = f0/sum(f0);
N_0 = N_T0*f0;                                                             

Xi  = xi_SR*exp(k_xi*(1 - R));
Xi  = Xi - diag(diag(Xi));

AA_aux = Xi' - diag(sum(Xi, 2));
b      = b_S*b_R./(b_R + r.^alpha_b*(b_S - b_R));
d_max  = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);

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
    
    N_T                    = sum(xout, 2);
    N_Tmod(1:ntexp, iexp) = N_T(texp_ind);

end

% ----------------------------------------------------------------------- %
% Value of cost function:
N_Tmod  = reshape(N_Tmod, [], 1);
N_Tave_data = reshape(N_Tave_data, [], 1);

diff_NT = N_Tave_data - N_Tmod;

J       = ntexp*Nexp*log(sum(diff_NT.^2./(N_Tmod.^var_b))/(ntexp*Nexp)) + var_b*sum(log(N_Tmod));

end
