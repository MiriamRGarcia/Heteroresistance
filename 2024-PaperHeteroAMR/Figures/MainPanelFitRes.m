%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MainPanelFitRes: Plot results of model calibration using MLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

addpath('Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of trajectories in the calibration problem:
m_traj = 3;

% Method to generate BD trajectories:
method = 'SSA';                                                            % 'SSA'; 'RSSA';

% Low detection limit:
LDL    = 10;

% Set ODE solver precision:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

% Save figure (=1) or not (=0):
fig_print = 0;

% Format to save figure:
fig_form  = '-dpng'; %'djpg';

% Define colors for plot (one per experiment):
cc1 = [39, 183, 222]/256;
cc2 = [76,57,87]/256;
cc3 = [129,23,27]/256;                                            
cc4 = [237,106,90]/256;
cc5 = [51,115,87]/256;

cc  = [cc1;cc2;cc3;cc4;cc5];

% Define markers on sampling times (one per experiment):
mks{1} = 's';                                                            
mks{2} = 'diamond';
mks{3} = 'v';
mks{4} = '^';
mks{5} = '>';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user-defined settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load common experimental setup:
file_name = sprintf('../PE/Results/resPE_MNHo_%utraj.mat', m_traj);
load(file_name, 'tsim', 'texp', 'r', 'Cexp')

m_r = numel(r);
m_t = numel(tsim);
m_e = numel(Cexp);

texp_ind = find(ismember(tsim, texp));
ntexp    = numel(texp_ind);


% ----------------------------------------------------------------------- %
% First subplot (MNHo case):
file_name = sprintf('../PE/Results/resPE_MNHo_%utraj.mat', m_traj);

load(file_name, 'pars', 'seed', 'sd')

% Generate random trajectories with exact parameter values:
rng(seed)

b_S       = pars(1);                                                       % Obtain exact parameter values;
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

f_0 = exp(-lambda_T0*r);                                                   % Initial condition;
f_0 = f_0/sum(f_0);
N_0 = N_T0*f_0;


N_T = zeros(m_t, m_e);


R  = repmat(r, 1, m_r) - repmat(r.', m_r, 1);                                % Modification rates matrix;                 
R  = R - triu(R) + tril(R).';   
Xi = xi_SR*exp(k_xi*(1 - R));
Xi = Xi - diag(diag(Xi));


AA_aux = Xi' - diag(sum(Xi, 2));                                           % Auxiliary coefficient matrix;


b     = b_S*b_R./(b_R + r.^alpha_b*(b_S - b_R));                           % Birth and maximal kill rates;
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);

for iexp = 1:m_e
    
    %-------------------------------------------------%
    % Calculate coefficient matrix of the state system:
    
    % Drug concentration:
    CC = Cexp(iexp);
    
    % Calculate kill rate at current time:
    HC = CC^H_d/(CC^H_d + EC_50d^H_d);
    d  = d_max*HC;
    
    % Coefficient matrix:
    AA = AA_aux + diag(b - d);
    
    %-------------------------------------------------%
    % Calculate population size with exact parameters:
    [~, xout] = ode15s(@(t,s) Odes_cte(t, s, AA), tsim, N_0, ODEoptions);
   
    % Total population:
    N_T(1:m_t, iexp) = sum(xout, 2);
 
end

% Trajectories with homoscedastic normal noise:
N_Tdata = zeros(m_t, m_e, m_traj);

for itraj = 1:m_traj
    N_Tdata(1:m_t, 1:m_e, itraj) = log10(N_T(1:m_t, 1:m_e)) + sd*randn(m_t, m_e);
end

% Load fit results:
load(file_name, 'pars_opt', 'N_Tave_data')

% Obtain optimal parameter values to plot model total counts:
b_S       = pars_opt(1);
b_R       = pars_opt(2);
alpha_b   = pars_opt(3);
d_maxS    = pars_opt(4);
beta_d    = pars_opt(5);
alpha_d   = pars_opt(6);
EC_50d    = pars_opt(7);
H_d       = pars_opt(8);
xi_SR     = pars_opt(9);
k_xi      = pars_opt(10);
N_T0      = pars_opt(11);
lambda_T0 = pars_opt(12);

% Initial condition:
f_0 = exp(-lambda_T0*r);
f_0 = f_0/sum(f_0);
N_0 = N_T0*f_0;

% Mutation rates matrix:
Xi = xi_SR*exp(k_xi*(1 - R));
Xi = Xi - diag(diag(Xi));

% Auxiliary coefficient matrix:
AA_aux = Xi' - diag(sum(Xi, 2));

% Calculate growth rate:
b = b_S*b_R./(b_R + r.^alpha_b*(b_S - b_R));

% Calculate death rate:
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);

% Initialice total population size:
N_Tmod = zeros(m_t, m_e);

% Initialise legend:
lgd = cell(m_e, 1);

% Figure:
fig = figure;

set(gcf,'color','w');

transp    = [0.3 0.5 0.7];
ind_tplot = 1:1000:m_t;

N_Tdata(N_Tdata < log10(LDL)) = log10(LDL);

% First subplot:
subplot(1, 3, 1)

hold on

for iexp = 1:m_e
    
    %-------------------------------------------------%
    % Calculate coefficient matrix of the state system:
    
    % Drug concentration:
    CC = Cexp(iexp);
    
    % Calculate kill rate at current time:
    HC     = CC^H_d/(CC^H_d + EC_50d^H_d);
    d    = d_max*HC;
    
    % Coefficient matrix:
    AA     = AA_aux + diag(b - d);
    
    %-------------------------------------------------%
    % Calculate population size:
    [~, xout] = ode15s(@(t,s) Odes_cte(t, s, AA), tsim, N_0, ODEoptions);
   
    % Total population:
    N_Tmod(1:m_t, iexp) = log10(sum(xout, 2));
    
    
    for itraj = 1:m_traj
        plot(tsim(ind_tplot), N_Tdata(ind_tplot, iexp, itraj), 'Color', [cc(iexp,:) transp(itraj)], 'HandleVisibility', 'off')
    end
    
    plot(tsim, N_Tmod(1:m_t, iexp), 'Color', cc(iexp, :),'LineWidth', 2)
    
    plot(texp, N_Tave_data(1:ntexp, iexp), mks{iexp}, 'Color', cc(iexp,:),'MarkerSize', 10, 'MarkerFaceColor', cc(iexp,:), 'HandleVisibility', 'off')
    
    lgd{iexp} = sprintf('$C=%.1f$ mg/L', Cexp(iexp));

end

plot(tsim, log10(LDL)*ones(m_t, 1), '--', 'Color', 'k', 'LineWidth', 1.5)
text(11, log10(LDL) - 0.5, 'LDL', 'Interpreter', 'Latex', 'FontSize', 17)

set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'Latex', 'XTick', [0 12 24 36 48])

ylabel({'$N_T$ $\log_{10}$(CFUS/mL)'}, 'Interpreter', 'Latex', 'FontSize', 17)

xlim([0 48])
ylim([0 20])

legend(lgd, 'Orientation', 'Horizontal', 'Interpreter', 'Latex', 'Color', 'w', 'Edgecolor', 'None', ...
    'FontSize', 17, 'Position',[0.0696557245796627 -0.0266609726788979 0.907926341072859 0.0795710255889508])

   
hold off
box off


% ----------------------------------------------------------------------- %
% Second subplot (MNHe case):
file_name = sprintf('../PE/Results/resPE_MNHe_%utraj.mat', m_traj);

load(file_name, 'pars_opt', 'pars_var', 'N_Tave_data', 'seed')

N_Tave_data = log10(N_Tave_data);

% Generate random trajectories:
rng(seed)

var_a = pars_var(1);
var_b = pars_var(2);

N_Tdata = zeros(m_t, m_e, m_traj);

% Heterocedastic variance of trajectories:
vartraj      = var_a*N_T(1:m_t, 1:m_e).^var_b;

 for itraj = 1:m_traj    
     aux = N_T(1:m_t, 1:m_e) + sqrt(vartraj).*randn(m_t, m_e);
     
     aux(aux < LDL) = LDL;
     
     N_Tdata(1:m_t, 1:m_e, itraj) = log10(aux);
 end
 
% Obtain optimal parameter values to plot model fit:
b_S       = pars_opt(1);
b_R       = pars_opt(2);
alpha_b   = pars_opt(3);
d_maxS    = pars_opt(4);
beta_d    = pars_opt(5);
alpha_d   = pars_opt(6);
EC_50d    = pars_opt(7);
H_d       = pars_opt(8);
xi_SR     = pars_opt(9);
k_xi      = pars_opt(10);
N_T0      = pars_opt(11);
lambda_T0 = pars_opt(12);

% Initial condition:
f_0 = exp(-lambda_T0*r);
f_0 = f_0/sum(f_0);
N_0 = N_T0*f_0;

% Mutation rates matrix:
Xi = xi_SR*exp(k_xi*(1 - R));
Xi = Xi - diag(diag(Xi));

% Auxiliary coefficient matrix:
AA_aux = Xi.' - diag(sum(Xi, 2));

% Calculate growth rate:
b = b_S*b_R./(b_R + r.^alpha_b*(b_S - b_R));

% Calculate death rate:
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);

% Initialice total population size:
N_Tmod = zeros(m_t, m_e);

ind_tplot = 1:1000:m_t;

subplot(1, 3, 2)

hold on
    
for iexp = 1:m_e
    
    %-------------------------------------------------%
    % Calculate coefficient matrix of the state system:
    
    % Drug concentration:
    CC = Cexp(iexp);
    
    % Calculate kill rate at current time:
    HC     = CC^H_d/(CC^H_d + EC_50d^H_d);
    d    = d_max*HC;
    
    % Coefficient matrix:
    AA     = AA_aux + diag(b - d);
    
    %-------------------------------------------------%
    % Calculate population size:
    [~, xout] = ode15s(@(t,s) Odes_cte(t, s, AA), tsim, N_0, ODEoptions);
    
    % Total population:
    N_Tmod(1:m_t, iexp) = log10(sum(xout, 2));
    
    for itraj = 1:m_traj
        plot(tsim(ind_tplot), N_Tdata(ind_tplot, iexp, itraj), 'Color', [cc(iexp,:) transp(itraj)], 'HandleVisibility', 'off')
    end
    
     plot(tsim, N_Tmod(1:m_t,iexp), 'Color', cc(iexp, :),'LineWidth', 2)
     plot(texp, N_Tave_data(1:ntexp, iexp), mks{iexp}, 'Color', cc(iexp,:),'MarkerSize', 10, 'MarkerFaceColor', cc(iexp,:), 'HandleVisibility', 'off')

end

plot(tsim, log10(LDL)*ones(m_t, 1), '--', 'Color', 'k', 'LineWidth', 1.5)
text(11, log10(LDL) - 0.5, 'LDL', 'Interpreter', 'Latex', 'FontSize', 17)

xlim([0 48])
ylim([0 20])


set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'Latex', 'XTick', [0 12 24 36 48])

xlabel({'','Time (h)'}, 'Interpreter', 'Latex', 'FontSize', 17)

hold off
box off

% ----------------------------------------------------------------------- %
% Third subplot (PN case):
file_name = sprintf('../PE/Results/resPE_PN_%utraj.mat', m_traj);

load(file_name, 'pars_opt', 'N_Tave_data', 'itraj')

count = 1;
for ii = itraj
    file_name = sprintf('../SSA/Results/res%s_%03u', method, ii);
    load(file_name, 'N_T')
    
    N_Tdata(1:m_t, 1:m_e, count) = N_T(1:m_t, 1:m_e);    
    count = count + 1;
end

N_Tdata(N_Tdata < LDL) = LDL;
N_Tdata     = log10(N_Tdata);
N_Tave_data = log10(N_Tave_data);

% Obtain optimal parameter values to plot model average:
b_S       = pars_opt(1);
b_R       = pars_opt(2);
alpha_b   = pars_opt(3);
d_maxS    = pars_opt(4);
beta_d    = pars_opt(5);
alpha_d   = pars_opt(6);
EC_50d    = pars_opt(7);
H_d       = pars_opt(8);
xi_SR     = pars_opt(9);
k_xi      = pars_opt(10);
N_T0      = pars_opt(11);
lambda_T0 = pars_opt(12);

% Initial condition:
f_0  = exp(-lambda_T0*r);
f_0  = f_0/sum(f_0);
N_0  = N_T0*f_0;

% Initialice total population size:
N_Tmod = zeros(m_t, m_e);

% Mutation rates matrix:
Xi = xi_SR*exp(k_xi*(1 - R));
Xi = Xi - diag(diag(Xi));

AA_aux = Xi' - diag(sum(Xi, 2));

% Calculate growth rate:
b = b_S*b_R./(b_R + r.^alpha_b*(b_S - b_R));

% Calculate death rate:
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);

% ----------------------------------------------------------------------- %
% Plot fit results with process noise:

ind_tplot = 1:m_t;

transp = [0.6 0.6 0.6];

subplot(1, 3, 3)

hold on

for iexp = 1:m_e
    
    %-------------------------------------------------%
    % Calculate coefficient matrix of the state system:
    
    % Drug concentration:
    CC = Cexp(iexp);
    
    % Calculate kill rate at current time:
    HC     = CC^H_d/(CC^H_d + EC_50d^H_d);
    d    = d_max*HC;
    
    % Coefficient matrix:
    AA     = AA_aux + diag(b - d);
    
    %-------------------------------------------------%
    % Calculate population size:
    [~, xout] = ode15s(@(t,s) Odes_cte(t, s, AA), tsim, N_0, ODEoptions);
   
    % Total population:
    N_Tmod(1:m_t, iexp) = log10(sum(xout, 2));
    
    for itraj = 1:m_traj
        plot(tsim(ind_tplot), N_Tdata(ind_tplot, iexp, itraj), 'Color', [cc(iexp,:) transp(itraj)], 'HandleVisibility', 'off')
    end
    
    aux = reshape(N_Tdata(1:m_t, iexp, 1), m_t, 1);
    
    NaN_ind = find(isnan(aux));
    
    N_Tmod(NaN_ind, iexp) = NaN;
    
    plot(tsim, N_Tmod(1:m_t, iexp), 'Color', cc(iexp, :),'LineWidth', 2)
    
    plot(texp, N_Tave_data(1:ntexp, iexp), mks{iexp}, 'Color', cc(iexp,:),'MarkerSize', 10, 'MarkerFaceColor', cc(iexp,:), 'HandleVisibility', 'off')

end
plot(tsim, log10(LDL)*ones(m_t, 1), '--', 'Color', 'k', 'LineWidth', 1.5)
text(11, log10(LDL) - 0.5, 'LDL', 'Interpreter', 'Latex', 'FontSize', 17)

xlim([0 48])
ylim([0 20])


set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'Latex', 'XTick', [0 12 24 36 48])
box off

% ----------------------------------------------------------------------- %
% Zoom into figure:

% Plot rectangle around zoom area:
x0z = 3.25;
xfz = 21.5;

y0z = 2.25;
yfz = 4.5;

npz = 100;
xx  = linspace(x0z - 0.3, xfz + 0.2, npz);
yy  = linspace(y0z - 0.2, yfz + 0.3, npz);

plot(xx, y0z*ones(npz,1) - 0.2, xx, yfz*ones(npz,1) + 0.3, 'LineWidth', 1.5, 'Color', 'k', 'HandleVisibility', 'off')
plot(x0z*ones(npz,1) - 0.3, yy, xfz*ones(npz,1) + 0.2, yy, 'LineWidth', 1.5, 'Color', 'k', 'HandleVisibility', 'off')

plot(linspace(xx(end), 46, npz), 0.5*(yy(1) + yy(end))*ones(1,npz), 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'k', 'HandleVisibility', 'off')
plot(46*ones(1, npz), linspace(0.5*(yy(1) + yy(end)), 9, npz), 'LineWidth', 1.5,  'LineStyle', '--', 'Color', 'k', 'HandleVisibility', 'off')
% create a new pair of axes inside current figure
%axes('FontSize',6, 'LineWidth', 5, 'position',[0.825520833333334,0.449764544456641,0.148437499999999,0.44825],...
%     'outerposition',[0.800621639784947,0.389264544456641,0.191532258064515,0.55])

hAx = axes('Position',[0.825520833333334 0.466051189407781 0.148437499999999 0.44825]);


box on % put box around new pair of axes

indexOfInterest = find(tsim < xfz & tsim > x0z); % range of t near perturbation

aux_ind  = intersect(indexOfInterest, texp_ind);
aux_ind1 = find(texp_ind > indexOfInterest(1), 1 );
aux_ind2 = find(texp_ind < indexOfInterest(end), 1, 'last');


for iexp = 1:m_e
    for itraj = 1:m_traj
        plot(tsim(indexOfInterest),N_Tdata(indexOfInterest, iexp, itraj), 'Color', [cc(iexp,:) transp(itraj)], 'HandleVisibility', 'off') % plot on new axes
        hold on
        plot(tsim(indexOfInterest),N_Tmod(indexOfInterest, iexp),'Color', cc(iexp,:),'LineWidth', 2, 'HandleVisibility', 'off') % plot on new axes
        
    end
    if numel(aux_ind1) > 0
        plot(tsim(aux_ind), N_Tave_data(aux_ind1:aux_ind2, iexp), mks{iexp}, 'Color', cc(iexp,:),'MarkerSize', 10, 'MarkerFaceColor', cc(iexp,:), 'HandleVisibility', 'off')
    end
end
axis tight
xlim([x0z xfz])
ylim([y0z yfz])

set(hAx,'TickLabelInterpreter','latex', 'YTick', y0z:0.3:yfz);
hAx.YAxis.TickLabelFormat = '%.1f';
hAx.LineWidth = 1;            % set the axis linewidth for box/ticks
hAx.FontSize  = 13;
%hAx.TickLabelFormat = '%.1f';

% ----------------------------------------------------------------------- %
% Save figure:
set(gcf, 'Position', [136         137        1590         823])

% ----------------------------------------------------------------------- %
% Print figure if desired:
if fig_print == 1
    print('-r720','PanelFitRes', fig_form)
end
