%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results of model calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotPE(tmod, r, par_opt, logNT_ave_data, col, mks) 

% Obtain parameter values:
bS        = par_opt(1);
bR        = par_opt(2);
alpha_b   = par_opt(3);
d_maxS    = par_opt(4);
alpha_d   = par_opt(5);
beta_d    = par_opt(6);
EC_50d    = par_opt(7);
H_d       = par_opt(8);
xi_SR     = par_opt(9);
k_xi      = par_opt(10);
N_T0      = par_opt(11);
lambda_T0 = par_opt(12);

% Initial condition:
f0  = exp(-lambda_T0*r);
f0  = f0/sum(f0);
N0  = N_T0*f0;

% Initialice coefficient matrix:
Xi = xi_SR*exp(k_xi*(1 - RR));
Xi = Xi - diag(diag(Xi));

AA_aux = Xi' - diag(sum(Xi, 2));

% Birth rate:
b     = bS*bR./(bR + r.^alpha_b*(bS - bR));

% Maximal death rate:
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);

% Plot results:
figure

hold on

lgd = cell(Nexp, 1);

for iexp = 1:Nexp
    
    C  = Cexp(iexp);
    HC = C^H_d/(C^H_d + EC_50d^H_d);
    d  = d_max*HC;

    AA = AA_aux + diag(b - d);

    % Solve ODEs system:
    [~, xout] = ode15s(@(s,y) Odes_cte(s, y, AA), tmod, N0, ODEoptions);
    
    NT_mod    = sum(xout, 2);
    logNT_mod = log10(NT_mod);
    
    plot(texp, logNT_ave_data, mks{iexp}, 'Color', col(iexp,:), 'MarkerSize', 10, 'MarkerFaceColor', col(iexp,:), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off')
    plot(tmod, logNT_mod, 'Color', col(iexp, :), 'LineWidth', 1.5)
    
    lgd{iexp} = sprintf('$C=%0.2f$ (mg/L)', C);

end

xlabel('Time (h)', 'Interpreter', 'Latex', 'FontSize', 10)
ylabel('Total Bacterial Count ($\log_{10}$CFUS/mL)', 'Interpreter', 'Latex', 'FontSize', 10)

legend(lgd, 'Location', 'Best', 'Interpreter', 'Latex')

hold off

box off

end
