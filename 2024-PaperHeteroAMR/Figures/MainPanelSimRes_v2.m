%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot panel with simulation results of the TKA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
addpath('Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Load data of Gillespie trajectories:

% Number of trajectories to load:
m_traj = 1000;

% Choose implementation (direct method = SSA or rejection based = RSSA):
method = 'RSSA'; % = 'SSA'; = 'RSSA';

% Load experimental setup:
file_name = sprintf('../SSA/Results/res%s_001.mat', method);

load(file_name, 'r', 'tsim', 'Cexp', 'par')

% Problem sizes:
m_r = numel(r);
m_t = numel(tsim);
m_e = numel(Cexp); 

% Initialise Gillespie trajectories of total cell counts (CFUS/mL):
CFUST_data = zeros(m_t, m_e, m_traj);

% Loop in the trajectories:
for itraj = 1:m_traj

    % Create file name:
    file_name = sprintf('../SSA/Results/res%s_%03u.mat', method, itraj);

    % Load results file:
    load(file_name, 'CFUS_exp', 'CFUST_exp')

    % Preprocesing data:
    %CFUST_exp(CFUST_exp == N_TL) = NaN;
    
    for iexp = 1:m_e
        aux       = [0;
                     CFUST_exp(1:(m_t - 1), iexp)];
                 
        [rind, ~] = find((CFUST_exp(1:m_t, iexp) - aux == 0) & (CFUST_exp(1:m_t, iexp) - repmat(CFUST_exp(end, iexp), m_t, 1) == 0) & (CFUST_exp(1:m_t, iexp) - repmat(max(CFUST_exp(1:m_t, iexp)), m_t, 1) == 0));
        
        CFUST_exp(rind, iexp)      = NaN;
        
        CFUS_exp(rind, 1:m_r, iexp) = NaN;
    end

    % Almacenate total cell count data:  
    CFUST_data(1:m_t, 1:m_e, itraj) = CFUST_exp;
 
end

% Sample average of total cell counts:
CFUST_ave_data = sum(CFUST_data, 3)/m_traj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Obtain median and inter-quartile range of data:
CFUST_data_median   = median(CFUST_data, 3);
CFUST_data_quarlow  = quantile(CFUST_data, 0.25, 3);
CFUST_data_quarup   = quantile(CFUST_data, 0.75, 3);
IQR                 = CFUST_data_quarup - CFUST_data_quarlow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Calculate deterministic average of the BD model:

% Set ODE solver precision:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

% Obtain parameter values:
bS        = par(1);
bR        = par(2);
alpha_b   = par(3);

d_maxS    = par(4);
beta_d    = par(5);
alpha_d   = par(6);

EC_50d    = par(7);
H_d       = par(8);

xi_SR     = par(9);
k_xi      = par(10);
    
N_T0      = par(11);
lambda_T0 = par(12);

% Initial condition:
f0  = exp(-lambda_T0*r);
f0  = f0/sum(f0);
N0  = N_T0*f0;

% Initialice population sizes:
CFUST_mod = zeros(m_t, m_e);

% Initialise coefficient matrix:
RR = repmat(r, 1, m_r) - repmat(r.', m_r, 1);                      
RR = RR - triu(RR) + tril(RR).';

Xi = xi_SR*exp(k_xi*(1 - RR));
Xi = Xi - diag(diag(Xi));

AA_aux = Xi.' - diag(sum(Xi, 2));

% Calculate growth rate at current time:
b = bS*bR./(bR + r.^alpha_b*(bS - bR));

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
    [~, xout] = ode15s(@(t, s) Odes_cte(t, s, AA), tsim, N0, ODEoptions);
    
    %-------------------------------------------------%
    % Almacenate cell numbers for current antimicrobial concentration:
    CFUST_mod(1:m_t, iexp)      = sum(xout, 2);

end

rmpath('Functions')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (4) Plot results (fitness cost and kill rate):
close all
addpath('Functions')

load(file_name, 'Cexp')

MIC_S = EC_50d*(bS/(d_maxS - bS))^(1/H_d);

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

fig = figure;

set(gcf, 'Position', [180 122 1544 854], 'Color', 'w')

nr_plot = 1e3;
rr_plot = linspace(0, 1, nr_plot);

% Calculate growth rate:
b_plot = bS*bR./(bR + rr_plot.^alpha_b*(bS - bR));

% ----------------------------------------------------------------------- %
% Plot birth-rate:
ss = subplot(2, 3, 1);
hold on
plot(rr_plot, b_plot, 'Color', bb5, 'LineWidth', 1.5)
plot(r, bS*ones(m_r, 1), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
plot(r, bR*ones(m_r, 1), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
text(0.1, bS + 0.02, '$b_S$', 'Interpreter', 'Latex', 'FontSize', 17)
text(0.1, bR + 0.02, '$b_R$', 'Interpreter', 'Latex', 'FontSize', 17)

rr_max = (bR*(alpha_b - 1)/((alpha_b + 1)*(bS - bR)))^(1/alpha_b);
plot(rr_max*ones(m_r, 1), linspace(bR, bS, m_r), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
ttb = text(rr_max, 0.32, '$r_b$', 'Interpreter', 'Latex', 'FontSize', 17);
ttb.Position = [0.661498708010336 0.333680781758958 0];

hold off

xlabel('$r$', 'Interpreter', 'Latex', 'FontSize', 17)
ylabel('$b(r)$ (1/h)', 'Interpreter', 'Latex', 'FontSize', 17)

xticks = [0 0.2 0.4 0.6 0.8 1.0];
set(ss, 'XTick', xticks, 'TickLabelInterpreter', 'Latex', 'FontSize', 14)
ss.YAxis.TickLabelFormat = '%,.2f';

% ----------------------------------------------------------------------- %
% Plot death rate:
Exps = [1 2 3 4 5 9];

dmax_plot = d_maxS*beta_d^alpha_d*(1 - rr_plot.^alpha_d)./(beta_d^alpha_d + rr_plot.^alpha_d);

ss = subplot(2, 3, 4);

hold on

for iexp = Exps
    CC     = Cexp(iexp);
    HC     = CC^H_d/(CC^H_d + EC_50d^H_d);
    d_plot = dmax_plot.*HC;
    plot(rr_plot, d_plot, 'Color', cc(iexp, :), 'LineWidth', 1.5)
end
plot(r, d_maxS*ones(m_r, 1), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
text(0.1, d_maxS + 0.3, '$d_{maxS}$', 'Interpreter', 'Latex', 'FontSize', 17)

rr_max = ((alpha_d - 1)*beta_d^alpha_d/(alpha_d + 1))^(1/alpha_d);
plot(rr_max*ones(m_r, 1), linspace(0, d_maxS, m_r), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
ttk = text(rr_max, -0.1, '$r_d$', 'Interpreter', 'Latex', 'FontSize', 17);
ttk.Position = [0.286472458455655 -0.217263843648209 0];
hold off

xlabel('$r$', 'Interpreter', 'Latex', 'FontSize', 17)
ylabel('$d(r,C)$ (1/h)', 'Interpreter', 'Latex', 'FontSize', 17)

set(ss, 'XTick', xticks, 'TickLabelInterpreter', 'Latex','FontSize', 14)
ss.YAxis.TickLabelFormat = '%,.2f';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5) Plot results (statistics: mean, interquartile range etc):

% Tranform data to log-scale:
logCFUST_mod      = log10(CFUST_mod);
logCFUST_data     = log10(CFUST_data);
logCFUST_ave_data = log10(CFUST_ave_data);
logCFUST_data_median  = log10(CFUST_data_median);
logCFUST_data_quarlow = log10(CFUST_data_quarlow);
logCFUST_data_quarup  = log10(CFUST_data_quarup);

ss = subplot(2, 3, [2 5]);

hold on

Exps     = [1 2 3 4 5 9];
Nexp_aux = numel(Exps);

lgd       = cell(Nexp_aux, 1);
lgd_count = 1;
MIC_label = [0 1 2 4 8 10 12 16 32];

for iexp = Exps
    
    NaNind = find(isnan(logCFUST_ave_data(1:m_t, iexp)));
    
    % ----------------------------------------------- %
    % Plot model average of data:
    logCFUST_mod(NaNind, iexp) = NaN;
    
    plot(tsim, logCFUST_mod(1:m_t, iexp), 'Color', cc(iexp, :), 'LineWidth', 1.0)
    
    % ----------------------------------------------- %
    % Plot average value of trajectories:
    plot(tsim, logCFUST_ave_data(1:m_t, iexp), 'LineStyle', '--', 'Color', cc(iexp, :), 'LineWidth', 1.0, 'HandleVisibility', 'off')

    % ----------------------------------------------- %
    % Shade area between quartile 0.25 and 0.75:
    lowlim = logCFUST_data_quarlow(1:m_t, iexp);
    uplim  = logCFUST_data_quarup(1:m_t, iexp);
    
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
    CFUST_data_aux = reshape(CFUST_data(1:m_t, iexp, 1:m_traj), m_t, m_traj);
    IQR_aux        = IQR(1:m_t, iexp);
    limlow         = repmat(CFUST_data_quarlow(1:m_t, iexp) - 1.5*IQR_aux, 1, m_traj);
    limup          = repmat(CFUST_data_quarup(1:m_t, iexp)  + 1.5*IQR_aux, 1, m_traj);
    outind         = find(limup < CFUST_data_aux | CFUST_data_aux < limlow);
     
    % ----------------------------------------------- %
    % Scatter plot with data outliers:
    %taux2 = repmat(tsim.', 1, Ntraj);
    
    %outind_plot = outind(1:10000:end);
    %scatter(taux2(outind_plot), log10(CFUST_data_aux(outind_plot)), [], cc(iexp,:), 'Marker', 'x', 'HandleVisibility', 'off')
    
    CFUST_data_aux(outind) = NaN;
    
    minCFUST_data = min(CFUST_data_aux, [], 2);
    maxCFUST_data = max(CFUST_data_aux, [], 2);

    if numel(NaNind) > 0
        minCFUST_data(NaNind:end) = [];
        maxCFUST_data(NaNind:end) = [];
    end
    
    lowlim = log10(minCFUST_data);
    uplim  = log10(maxCFUST_data);
    
    coord_low = [taux.', lowlim];
    coord_up  = [taux.', uplim];

    coord_combine = [coord_up;flipud(coord_low)];
    
    fill(coord_combine(:,1),coord_combine(:,2), cc(iexp, :), 'EdgeColor', 'None', 'FaceAlpha', 0.2, 'HandleVisibility', 'off')
    
  
    % Build legend:
    lgd{lgd_count} = sprintf('$C=%u MIC_S$', MIC_label(iexp));%strcat('$C=\;$ ', num2str(Cexp(iexp)), ' (mg/L)');
    
    lgd_count      = lgd_count + 1;

end

% Low detection limit:
%LDL = log10(10);

%plot(tsim, LDL*ones(nt,1), 'k', 'LineStyle', '--', 'LineWidth', 1.5)
xlabel('Time (h)', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$N_T$ ($log_{10}$CFUS/mL)', 'Interpreter', 'Latex', 'FontSize', 15)

xlim([0 40])

set(gca, 'XTick', [0.0 5.0 10.0 15.0 20.0 25.0 30.0 35.0 40.0])

legend(lgd, 'Orientation', 'Horizontal', 'Interpreter', 'Latex', 'edgecolor', 'none', 'FontSize', 15, 'Position',[0.183997965068409 0.00912237826224061 0.676653187471205 0.0297365115607248])
%text(10, LDL + 0.1, 'LDL', 'Interpreter', 'Latex', 'FontSize', 17)

hold off
box off
set(ss, 'TickLabelInterpreter', 'Latex', 'FontSize', 14)
%set(ss, 'FontSize', 13)
%ss.XAxis.TickLabelFormat = '%,.1f';
%ss.YAxis.TickLabelFormat = '%,.1f';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (6) Plot extinction probability:
load(file_name, 'Cexp')
Cexp = linspace(Cexp(1), Cexp(end-1), 50);
m_e   = numel(Cexp);
Pext = zeros(m_t, m_e);

% Initial condition for ODEs:
z0 = [N0;
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
%lCB  = get(CB,'Limits');
%tCB  = linspace(lCB(1), lCB(2), 4);
    
CB.FontSize = 13;
    
%set(CB,'Ticks',tCB,'TickLabelInterpreter', 'Latex')
     
%TLCB = arrayfun(@(x) sprintf('%.2g',x), tCB, 'un', 0);
      
%set(CB,'TickLabels',TLCB)
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

print('-r1080', 'PanelSimResHD_v2', '-dpng')

rmpath('Functions')
