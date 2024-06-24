%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MainPanelSimRes: Main file to plot simulation results of TKA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

addpath('Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of trajectories in the calibration problem:
m_traj = 1000;

% Choose implementation (direct method = SSA or rejection based = RSSA):
method = 'SSA';                                                            % = 'SSA'; = 'RSSA';

% Set ODE solver precision:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

% Save figure (=1) or not (=0):
fig_print = 0;

% Format to save figure:
fig_form  = '-dpng'; %'djpg';

% Define colors (brown colormap):
bb1 = [92,47,24]/256;
bb3 = [196,147,88]/256;
bb2 = [184,97,52]/256;
bb4 = [124,42,11]/256;
bb5 = [132,114,68]/256;
bb6 = [215,176,48]/256;
bb7 = [215,176,48]/256;
bb8 = [215,176,48]/256;
bb9 = [215,176,48]/256;

% Define colors (contrasting colormap):
cc1 = [39, 183, 222]/256;
cc2 = [76,57,87]/256;
cc3 = [129,23,27]/256;                                            
cc4 = [237,106,90]/256;
cc5 = [51,115,87]/256;
cc6 = [250,163,0]/256;
cc7 = [250,163,0]/256;
cc9 = [250,163,0]/256;
cc8 = [250,163,0]/256;
cc10 = [250,163,0]/256;

cc = [cc1;cc2;cc3;cc4;cc5;cc6;cc7;cc8;cc9;cc10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user-defined settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% Load data of Gillespie trajectories:

file_name = '../SSA/Results/resSSA_001.mat';
load(file_name, 'r', 'tsim', 'Cexp', 'pars')

% Problem sizes:
m_r = numel(r);
m_t = numel(tsim);
m_e = numel(Cexp); 

% Initialise Gillespie trajectories of total cell counts (CFUS/mL):
N_Tdata = zeros(m_t, m_e, m_traj);

% Loop in the trajectories:
for itraj = 1:m_traj

    % Create file name:
    file_name = sprintf('../SSA/Results/res%s_%03u.mat', method, itraj);

    % Load results file:
    load(file_name, 'N', 'N_T')

    % Almacenate total cell count data:  
    N_Tdata(1:m_t, 1:m_e, itraj) = N_T;
 
end

% Sample average of total cell counts:
CFUST_ave_data = sum(N_Tdata, 3)/m_traj;

% ----------------------------------------------------------------------- %
% Obtain median and inter-quartile range of data:
N_Tdata_median   = median(N_Tdata, 3);
N_Tdata_quarlow  = quantile(N_Tdata, 0.25, 3);
N_Tdata_quarup   = quantile(N_Tdata, 0.75, 3);
IQR              = N_Tdata_quarup - N_Tdata_quarlow;

% ----------------------------------------------------------------------- %
% Calculate deterministic average of the BD model:

% Obtain parameter values:
b_S       = pars(1);
b_R       = pars(2);
alpha_b   = pars(3);

d_maxS    = pars(4);
beta_d    = pars(5);
alpha_d   = pars(6);

EC_50d    = pars(7);
H_d       = pars(8);

xi_SR     = pars(9);
k_xi      = pars(10);
    
N_T0      = pars(11);
lambda_T0 = pars(12);

% Initial condition:
f_0  = exp(-lambda_T0*r);
f_0  = f_0/sum(f_0);
N_0  = N_T0*f_0;

% Initialice population sizes:
N_Tmod = zeros(m_t, m_e);

% Initialise coefficient matrix:
R = repmat(r, 1, m_r) - repmat(r.', m_r, 1);                      
R = R - triu(R) + tril(R).';

Xi = xi_SR*exp(k_xi*(1 - R));
Xi = Xi - diag(diag(Xi));

AA_aux = Xi.' - diag(sum(Xi, 2));

% Calculate growth rate at current time:
b = b_S*b_R./(b_R + r.^alpha_b*(b_S - b_R));

% Calculate kill rate at current time:
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);

for iexp = 1:m_e
    
    %-------------------------------------------------%
    % Calculate coefficient matrix of the state system:
    
    % Drug concentration:
    CC = Cexp(iexp);
    
    % Calculate kill rate at current time:
    HC = CC^H_d/(CC^H_d + EC_50d^H_d);
    d  = d_max*HC;

    AA = AA_aux + diag(b - d);
    
    %-------------------------------------------------%
    % Call to ODEs:
    [~, xout] = ode15s(@(t, s) Odes_cte(t, s, AA), tsim, N_0, ODEoptions);
    
    %-------------------------------------------------%
    % Almacenate cell numbers for current antimicrobial concentration:
    N_Tmod(1:m_t, iexp)      = sum(xout, 2);

end

% ----------------------------------------------------------------------- %
% Plot results (fitness cost and kill rate):
MIC_S = EC_50d*(b_S/(d_maxS - b_S))^(1/H_d);

fig = figure;

set(gcf, 'Position', [180 122 1544 854], 'Color', 'w')

nr_plot = 1e3;
r_plot  = linspace(0, 1, nr_plot);

b_plot  = b_S*b_R./(b_R + r_plot.^alpha_b*(b_S - b_R));

% Plot birth-rate:
ss = subplot(2, 3, 1);
hold on
plot(r_plot, b_plot, 'Color', bb5, 'LineWidth', 1.5)
plot(r, b_S*ones(m_r, 1), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
plot(r, b_R*ones(m_r, 1), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
text(0.1, b_S + 0.02, '$b_S$', 'Interpreter', 'Latex', 'FontSize', 17)
text(0.1, b_R + 0.02, '$b_R$', 'Interpreter', 'Latex', 'FontSize', 17)

r_max = (b_R*(alpha_b - 1)/((alpha_b + 1)*(b_S - b_R)))^(1/alpha_b);
plot(r_max*ones(m_r, 1), linspace(b_R, b_S, m_r), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
ttb = text(r_max, 0.32, '$r_b$', 'Interpreter', 'Latex', 'FontSize', 17);
ttb.Position = [0.661498708010336 0.333680781758958 0];

hold off

xlabel('$r$', 'Interpreter', 'Latex', 'FontSize', 17)
ylabel('$b(r)$ (1/h)', 'Interpreter', 'Latex', 'FontSize', 17)

xticks = [0 0.2 0.4 0.6 0.8 1.0];
set(ss, 'XTick', xticks, 'TickLabelInterpreter', 'Latex', 'FontSize', 14)
ss.YAxis.TickLabelFormat = '%,.2f';

% Plot death rate:
Exps = [1 2 3 4 5 9];

dmax_plot = d_maxS*beta_d^alpha_d*(1 - r_plot.^alpha_d)./(beta_d^alpha_d + r_plot.^alpha_d);

ss = subplot(2, 3, 4);

hold on

for iexp = Exps
    CC     = Cexp(iexp);
    HC     = CC^H_d/(CC^H_d + EC_50d^H_d);
    d_plot = dmax_plot.*HC;
    plot(r_plot, d_plot, 'Color', cc(iexp, :), 'LineWidth', 1.5)
end
plot(r, d_maxS*ones(m_r, 1), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
text(0.1, d_maxS + 0.3, '$d_{maxS}$', 'Interpreter', 'Latex', 'FontSize', 17)

r_max = ((alpha_d - 1)*beta_d^alpha_d/(alpha_d + 1))^(1/alpha_d);
plot(r_max*ones(m_r, 1), linspace(0, d_maxS, m_r), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
ttk = text(r_max, -0.1, '$r_d$', 'Interpreter', 'Latex', 'FontSize', 17);
ttk.Position = [0.286472458455655 -0.217263843648209 0];
hold off

xlabel('$r$', 'Interpreter', 'Latex', 'FontSize', 17)
ylabel('$d(r,C)$ (1/h)', 'Interpreter', 'Latex', 'FontSize', 17)

set(ss, 'XTick', xticks, 'TickLabelInterpreter', 'Latex','FontSize', 14)
ss.YAxis.TickLabelFormat = '%,.2f';

% ----------------------------------------------------------------------- %
% Plot results (statistics: mean, interquartile range etc):

% Tranform data to log-scale:
logN_Tmod          = log10(N_Tmod);
logN_Tdata         = log10(N_Tdata);
logN_Tave_data     = log10(CFUST_ave_data);
logN_Tdata_median  = log10(N_Tdata_median);
logN_Tdata_quarlow = log10(N_Tdata_quarlow);
logN_Tdata_quarup  = log10(N_Tdata_quarup);

ss = subplot(2, 3, [2 5]);

hold on

Exps     = [1 2 3 4 5 9];
Nexp_aux = numel(Exps);

lgd       = cell(Nexp_aux, 1);
lgd_count = 1;
MIC_label = [0 1 2 4 8 10 12 16 32];

for iexp = Exps
    
    NaNind = find(isnan(logN_Tave_data(1:m_t, iexp)));
    
    % ----------------------------------------------- %
    % Plot model average of data:
    logN_Tmod(NaNind, iexp) = NaN;
    
    plot(tsim, logN_Tmod(1:m_t, iexp), 'Color', cc(iexp, :), 'LineWidth', 1.0)
    
    % ----------------------------------------------- %
    % Plot average value of trajectories:
    plot(tsim, logN_Tave_data(1:m_t, iexp), 'LineStyle', '--', 'Color', cc(iexp, :), 'LineWidth', 1.0, 'HandleVisibility', 'off')

    % ----------------------------------------------- %
    % Shade area between quartile 0.25 and 0.75:
    lowlim = logN_Tdata_quarlow(1:m_t, iexp);
    uplim  = logN_Tdata_quarup(1:m_t, iexp);
    
    if numel(NaNind) > 0
        lowlim(NaNind:end) = [];
        uplim(NaNind:end)  = [];
    
        taux = tsim(1:NaNind - 1);
    else
        taux = tsim;
    end
    
    coord_up  = [taux.', uplim];
    coord_low = [taux.', lowlim];
    
    coord_combine = [coord_up;flipud(coord_low)];
    
    fill(coord_combine(:,1),coord_combine(:,2), cc(iexp, :), 'EdgeColor', 'None', 'FaceAlpha', 0.4, 'HandleVisibility', 'off')
    

    % ----------------------------------------------- %
    % Maximum an minimum of data without outliers:
    N_Tdata_aux = reshape(N_Tdata(1:m_t, iexp, 1:m_traj), m_t, m_traj);
    IQR_aux     = IQR(1:m_t, iexp);
    limlow      = repmat(N_Tdata_quarlow(1:m_t, iexp) - 1.5*IQR_aux, 1, m_traj);
    limup       = repmat(N_Tdata_quarup(1:m_t, iexp)  + 1.5*IQR_aux, 1, m_traj);
    outind      = find(limup < N_Tdata_aux | N_Tdata_aux < limlow);

    N_Tdata_aux(outind) = NaN;
    
    minN_Tdata = min(N_Tdata_aux, [], 2);
    maxN_Tdata = max(N_Tdata_aux, [], 2);

    if numel(NaNind) > 0
        minN_Tdata(NaNind:end) = [];
        maxN_Tdata(NaNind:end) = [];
    end
    
    lowlim = log10(minN_Tdata);
    uplim  = log10(maxN_Tdata);
    
    coord_low = [taux.', lowlim];
    coord_up  = [taux.', uplim];

    coord_combine = [coord_up;flipud(coord_low)];
    
    fill(coord_combine(:,1),coord_combine(:,2), cc(iexp, :), 'EdgeColor', 'None', 'FaceAlpha', 0.2, 'HandleVisibility', 'off')
    
  
    % Build legend:
    lgd{lgd_count} = sprintf('$C=%u MIC_S$', MIC_label(iexp));%strcat('$C=\;$ ', num2str(Cexp(iexp)), ' (mg/L)');
    
    lgd_count      = lgd_count + 1;

end

xlabel('Time (h)', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$N_T$ ($log_{10}$CFUS/mL)', 'Interpreter', 'Latex', 'FontSize', 15)

xlim([0 40])

set(gca, 'XTick', [0.0 5.0 10.0 15.0 20.0 25.0 30.0 35.0 40.0])

legend(lgd, 'Orientation', 'Horizontal', 'Interpreter', 'Latex', 'edgecolor', 'none', 'FontSize', 15, 'Position',[0.183997965068409 0.00912237826224061 0.676653187471205 0.0297365115607248])

hold off
box off
set(ss, 'TickLabelInterpreter', 'Latex', 'FontSize', 14)


% ----------------------------------------------------------------------- %
% Plot extinction probability:
Cexp = linspace(Cexp(1), Cexp(end-1), 50);
m_e   = numel(Cexp);
Pext = zeros(m_t, m_e);

% Initial condition for ODEs:
z0 = [N_0;
      0;
      0];

for iC = 1:m_e
    
    %-------------------------------------------------%
    % Calculate coefficient matrix of the state system:
    
    % Drug concentration:
    CC = Cexp(iC);
    
    % Calculate kill rate at current time:
    d_max  = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);
    HC     = CC^H_d/(CC^H_d + EC_50d^H_d);
    d      = d_max*HC;

    AA     = AA_aux + diag(b - d);
    
    %-------------------------------------------------%
    % Call to ODEs with extinction probability:
    [~, zout] = ode15s(@(t,s) OdesPext(t, s, b, d, AA), tsim, z0, ODEoptions);
    zz        = zout(1:m_t, m_r + 2);
    
    % Extinction probability:
    Pext(1:m_t, iC) = (zz./(1 + zz)).^N_T0;
end

Cticks = 0:0.8:6.4;%[0 0.2 0.8 1.6 6.4];

subplot(1, 3, 3)

% Call to contour plot:
contourf(tsim, Cexp, Pext.', 50, 'edgecolor', 'none')

% Colorbar of subplot:
CB   = colorbar;

CB.FontSize = 13;

set(CB, 'TickLabelInterpreter', 'Latex')
 
% Figure settings:
xlabel('Time (h)', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$C$ (mg/L)', 'Interpreter','Latex','FontSize', 15);

xlim([tsim(1) 40])
ylim([Cexp(1) Cexp(end)])

% Axis properties:
set(gca,'XTick',[0 10 20 30 40 48], 'YTick', Cticks, 'TickLabelInterpreter', 'Latex', 'FontSize', 14)

% Set ticks on y axis:
ax = gca;
ax.YAxis.TickLabelFormat = '%.2f';
    
box off

if fig_print == 1
    print('-r720', 'PanelSimRes', fig_form)
end

rmpath('Functions')
