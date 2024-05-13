%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PlotPE: Plot results of model calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotPE(tmod, texp, r, Cexp, par_opt, logN_Tave_data, col, mks, ODEoptions) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% tmod           = Simulation times for solve ODEs (nt x 1);
% r              = Discretisation of the AMR level (values between entire 
%                  sensitivity S = 0 and entire resistance R = 1) (nr x 1),
% par_opt        = Calibrated values of model parameters;
% logN_Tave_data = Total count data (log10 scale);
% col            = Matrix with colors for each experiment;
% mks            = Cell array with markers for each experiment;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain parameter values:
b_S       = par_opt(1);
b_R       = par_opt(2);
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
N_0 = N_T0*f0;

% Initialice coefficient matrix:
nr = numel(r);
R  = repmat(r, 1, nr) - repmat(r.', nr, 1);
R  = R - triu(R) + tril(R).';
Xi = xi_SR*exp(k_xi*(1 - R));
Xi = Xi - diag(diag(Xi));

AA_aux = Xi' - diag(sum(Xi, 2));

% Birth rate:
b     = b_S*b_R./(b_R + r.^alpha_b*(b_S - b_R));

% Maximal death rate:
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);

% Plot results:
figure

hold on

Nexp = numel(Cexp);
lgd  = cell(Nexp, 1);

for iexp = 1:Nexp
    
    C  = Cexp(iexp);
    HC = C^H_d/(C^H_d + EC_50d^H_d);
    d  = d_max*HC;

    AA = AA_aux + diag(b - d);

    % Solve ODEs system:
    [~, xout] = ode15s(@(s,y) Odes_cte(s, y, AA), tmod, N_0, ODEoptions);
    
    N_Tmod    = sum(xout, 2);
    logN_Tmod = log10(N_Tmod);
    
    plot(texp, logN_Tave_data(:,iexp), mks{iexp}, 'Color', col(iexp,:), 'MarkerSize', 10, 'MarkerFaceColor', col(iexp,:), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off')
    plot(tmod, logN_Tmod, 'Color', col(iexp, :), 'LineWidth', 1.5)
    
    lgd{iexp} = sprintf('$C=%0.2f$ (mg/L)', C);

end

xlabel('Time (h)', 'Interpreter', 'Latex', 'FontSize', 10)
ylabel('Total Bacterial Count ($\log_{10}N_T$)', 'Interpreter', 'Latex', 'FontSize', 10)

legend(lgd, 'Location', 'Best', 'Interpreter', 'Latex')

hold off

box off

end
