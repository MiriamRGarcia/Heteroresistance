%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost function for MLE in MNHe case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = costFun_MNHe(decVars, r, R, tmod, texp_ind, Cexp, NT_data, pars_nom,Var_data, Weights, ODEoptions)

% ----------------------------------------------------------------------- %
% Obtain parameter values:
if size(pars_nom) == size(decVars)
    decVars = pars_nom.*decVars;
else
    decVars = pars_nom.'.*decVars;
end

bS        = decVars(1);
bR        = decVars(2);
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
N0  = N_T0*f0;                                                             

Xi  = xi_SR*exp(k_xi*(1 - R));
Xi  = Xi - diag(diag(Xi));

AA_aux = Xi' - diag(sum(Xi, 2));
b      = bS*bR./(bR + r.^alpha_b*(bS - bR));
d_max  = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);

% ----------------------------------------------------------------------- %
% Solve ODEs:
Nexp   = numel(Cexp);
ntexp  = numel(texp_ind);

NT_mod = zeros(ntexp, Nexp);

for iexp = 1:Nexp
    
    C  = Cexp(iexp);
    d  = d_max*C^H_d/(C^H_d + EC_50d^H_d);

    AA = AA_aux + diag(b - d);
    
    [~, xout] = ode15s(@(s,y) Odes_cte(s, y, AA), tmod, N0, ODEoptions);
    
    NT                    = sum(xout, 2);
    NT_mod(1:ntexp, iexp) = log10(NT(texp_ind));

end

% ----------------------------------------------------------------------- %
% Value of cost function:
NT_mod  = reshape(NT_mod, [], 1);
NT_data = reshape(NT_data, [], 1);

diff_NT = NT_data - NT_mod;

J       = ntexp*Nexp*log(sum(diff_NT.^2./(NT_mod.^var_b))/(ntexp*Nexp)) + var_b*sum(log(NT_mod));

end