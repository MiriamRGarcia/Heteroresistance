%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% costFun_MNHe: Cost function for MLE (PN case)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = costFun_PN(decVars, r, R, tsim, texp_ind, Cexp, N_Tave_data,...
                        pars_nom, Var_data, Weights, ODEoptions, sim_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% decVars     = Decision variables (parameters to calibrate) (size: m_p x 1);
% r           = Discretisation of the AMR level (values between entire 
%               sensitivity 0 and entire resistance 1) (size: m_r x 1),
% R           = Auxiliary matrix with jumps in AMR level (size: m_r x m_r)
% tsim        = Simulation times for solving ODEs (m_t x 1);
% texp_ind    = Indexes of sampling times within tsim (m_texp x 1);
% Cexp        = Array of antimicrobial concentrations (assumed  
%               constant at each experiment) (size: m_e x 1),
% N_Tave_data = Total count data at sampling times for MLE (average of
%               m_traj replicates) (size: m_texp x 1),
% pars_nom    = Nomimal values of parameters for scaling the problem;
% Var_data    = Variance in total count data at sampling times (for PN case);
% Weights     = Weights to remove NaN data (for PN case);
% ODEoptions  = Options for ODE solver;
% sim_name    = Name of the function simulating the average BD model;
%
% OUTPUT:
% J = Value of the cost function;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem sizes:
m_e = numel(Cexp);

% ----------------------------------------------------------------------- %
% Obtain parameter values (remove scaling):
if size(pars_nom, 2) - size(decVars, 2) < 0
    decVars = pars_nom.'.*decVars;
else
    decVars = pars_nom.*decVars;
end

% ----------------------------------------------------------------------- %
% Calculate total average counts:
[~, N_Tmod] = feval(sim_name, r, R, tsim, Cexp, decVars, ODEoptions);

N_Tmod = N_Tmod(texp_ind, 1:m_e);

% ----------------------------------------------------------------------- %
% Calculate cost function (weighted least squares):
J  = sum(sum(Weights.*((N_Tave_data - N_Tmod).^2)./Var_data));



end
