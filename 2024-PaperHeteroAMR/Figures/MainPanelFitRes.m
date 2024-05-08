%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot calibration results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

addpath('Functions')

% Set ODE solver precision:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) General experimental setup for model calibration:

% Number of trajectories in the calibration problem:
Ntraj = 3;

% Low detection limit:
LDL   = 10;

% Load experimental setup:
load('C:\Users\RYZEN9\Desktop\NEREA\Matlab\2023-HeteroRes\Paper_IdentObs\PE\TP_data\resFit_TP_3traj', 'tmod', 'rr', 'par')

% Number of experiments:
bS     = par(1);
d_maxS = par(4);
EC_50d    = par(7);
H_d    = par(8);

MIC_S  = EC_50d*(bS/(d_maxS - bS))^(1/H_d);
Cexp   = MIC_S*[0.0 1.0 2.0 4.0 8.0].';
Nexp   = numel(Cexp);

% Time discretisation:
nt       = numel(tmod);
t0       = tmod(1); 
tf       = tmod(nt);             
texp     = [2 4 6 8 10 12 16 20 24 36 48]; 
texp_ind = find(ismember(tmod, texp));   
ntexp    = numel(texp);

% Discretisation of AMR level:
nr       = numel(rr);
ra       = rr(1);
rb       = rr(nr);
RR       = repmat(rr, 1, nr) - repmat(rr.', nr, 1);                       
RR       = RR - triu(RR) + tril(RR).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) First subplot: results with measurement homocesdastic noise
path_name = 'C:\Users\RYZEN9\Desktop\NEREA\Matlab\2023-HeteroRes\Paper_IdentObs\PE\DT_data\Homo\Log';
file_name = strcat(path_name, filesep, 'resFit_TD_', num2str(Ntraj),'traj_log.mat');

load(file_name, 'par', 'seed')

% ----------------------------------------------------------------------- %
% Generate random trajectories with exact parameter values:
rng(seed)

% Obtain exact parameter values to generate data:
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

% Standar deviation of measurement error:
sd = 0.5;

% Initial condition:
f0  = exp(-lambda_T0*rr);
f0  = f0/sum(f0);
N0  = N_T0*f0;

% Initialice total population size with exact parameter values:
NT = zeros(nt, Nexp);

% Mutation rates matrix:
Xi = xi_SR*exp(k_xi*(1 - RR));
Xi = Xi - diag(diag(Xi));

% Auxiliary coefficient matrix:
AA_aux = Xi' - diag(sum(Xi, 2));

% Calculate growth rate:
b = bS*bR./(bR + rr.^alpha_b*(bS - bR));

% Calculate maximal kill rate:
d_max = d_maxS*beta_d^alpha_d*(1 - rr.^alpha_d)./(beta_d^alpha_d + rr.^alpha_d);

for iexp = 1:Nexp
    
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
    [~, xout] = ode15s(@(t,s) Odes_cte(t, s, AA), tmod, N0, ODEoptions);
   
    % Total population:
    NT(1:nt, iexp) = sum(xout, 2);
 
end

% Trajectories with homoscedastic normal noise:
logCFUS_traj = zeros(nt, Nexp, Ntraj);

for itraj = 1:Ntraj
    logCFUS_traj(1:nt, 1:Nexp, itraj) = log10(NT(1:nt, 1:Nexp)) + sd*randn(nt, Nexp);
end

% ----------------------------------------------------------------------- %
% Obtain fit results:
load(file_name, 'optPar', 'logCFUS_data')

% Obtain optimal parameter values to plot model total counts:
optpar  = optPar(end, :).';

bS        = optpar(1);
bR        = optpar(2);
alpha_b   = optpar(3);
d_maxS    = optpar(4);
beta_d    = optpar(5);
alpha_d   = optpar(6);
EC_50d    = optpar(7);
H_d       = optpar(8);
xi_SR     = optpar(9);
k_xi      = optpar(10);
N_T0      = optpar(11);
lambda_T0 = optpar(12);

% Initial condition:
f0  = exp(-lambda_T0*rr);
f0  = f0/sum(f0);
N0  = N_T0*f0;

% Mutation rates matrix:
Xi = xi_SR*exp(k_xi*(1 - RR));
Xi = Xi - diag(diag(Xi));

% Auxiliary coefficient matrix:
AA_aux = Xi' - diag(sum(Xi, 2));

% Calculate growth rate:
b = bS*bR./(bR + rr.^alpha_b*(bS - bR));

% Calculate death rate:
d_max = d_maxS*beta_d^alpha_d*(1 - rr.^alpha_d)./(beta_d^alpha_d + rr.^alpha_d);

% Initialice total population size:
logCFUS_mod = zeros(nt, Nexp);

% ----------------------------------------------------------------------- %
% Plot fit results:

% Define colors for plot (contrasting colormap):
cc1 = [39, 183, 222]/256;
cc2 = [76,57,87]/256;
cc3 = [129,23,27]/256;                                            
cc4 = [237,106,90]/256;
cc5 = [51,115,87]/256;

cc = [cc1;cc2;cc3;cc4;cc5];

% Define markers on sampling times:
mks{1} = 's';                                                            
mks{2} = 'diamond';
mks{3} = 'v';
mks{4} = '^';
mks{5} = '>';

% Initialise legend:
lgd = cell(Nexp, 1);

% Figure:
fig = figure;

set(gcf,'color','w');

transp    = [0.3 0.5 0.7];
ind_tplot = 1:1000:nt;

logCFUS_traj(logCFUS_traj < log10(LDL)) = log10(LDL);

% First subplot:
subplot(1, 3, 1)

hold on

for iexp = 1:Nexp
    
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
    [~, xout] = ode15s(@(t,s) Odes_cte(t, s, AA), tmod, N0, ODEoptions);
   
    % Total population:
    logCFUS_mod(1:nt, iexp) = log10(sum(xout, 2));
    
    
    for itraj = 1:Ntraj
        plot(tmod(ind_tplot), logCFUS_traj(ind_tplot, iexp, itraj), 'Color', [cc(iexp,:) transp(itraj)], 'HandleVisibility', 'off')
    end
    
    plot(tmod, logCFUS_mod(1:nt, iexp), 'Color', cc(iexp, :),'LineWidth', 2)
    
    plot(texp, logCFUS_data(1:ntexp, iexp), mks{iexp}, 'Color', cc(iexp,:),'MarkerSize', 10, 'MarkerFaceColor', cc(iexp,:), 'HandleVisibility', 'off')
    
    lgd{iexp} = sprintf('$C=%.1f$ mg/L', Cexp(iexp));

end

plot(tmod, log10(LDL)*ones(nt, 1), '--', 'Color', 'k', 'LineWidth', 1.5)
text(11, log10(LDL) - 0.5, 'LDL', 'Interpreter', 'Latex', 'FontSize', 17)

set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'Latex', 'XTick', [0 12 24 36 48])

ylabel({'Total counts $\log_{10}$(CFUS/mL)'}, 'Interpreter', 'Latex', 'FontSize', 17)

xlim([0 48])
ylim([0 20])

legend(lgd, 'Orientation', 'Horizontal', 'Interpreter', 'Latex', 'Color', 'w', 'Edgecolor', 'None', ...
    'FontSize', 17, 'Position',[0.0696557245796627 -0.0266609726788979 0.907926341072859 0.0795710255889508])

   
hold off
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Plot fit results with measurement heteroscedastic noise
clear logCFUS_data logCFUS_mod logCFUS_traj optpar

path_name = 'C:\Users\RYZEN9\Desktop\NEREA\Matlab\2023-HeteroRes\Paper_IdentObs\PE\DT_data\Hetero';
file_name = strcat(path_name, filesep, 'resFit_TD_Hetero_', num2str(Ntraj),'traj_b2.mat');

load(file_name, 'par', 'seed')

% Exact parameters of heterocedastic variance:
vara = 1;
varb = par(13);


% ----------------------------------------------------------------------- %
% Generate random trajectories:
rng(seed)

logCFUS_traj = zeros(nt, Nexp, Ntraj);

% Heterocedastic variance of trajectories:
vartraj      = vara*NT(1:nt, 1:Nexp).^varb;

 for itraj = 1:Ntraj    
     aux = NT(1:nt, 1:Nexp) + sqrt(vartraj).*randn(nt, Nexp);
     
     aux(aux < LDL) = LDL;
     
     logCFUS_traj(1:nt, 1:Nexp, itraj) = log10(aux);
 end
 
% ----------------------------------------------------------------------- %
% Plot results with heteroscedastic normal noise:

load(file_name, 'CFUS_data', 'CFUS_mod')

logCFUS_data = log10(CFUS_data);
logCFUS_mod  = log10(CFUS_mod);

ind_tplot = 1:1000:nt;

subplot(1, 3, 2)

hold on
    
for iexp = 1:Nexp

    for itraj = 1:Ntraj
        plot(tmod(ind_tplot), logCFUS_traj(ind_tplot, iexp, itraj), 'Color', [cc(iexp,:) transp(itraj)], 'HandleVisibility', 'off')
    end
    
     plot(tmod, logCFUS_mod(1:nt,iexp), 'Color', cc(iexp, :),'LineWidth', 2)
     plot(texp, logCFUS_data(1:ntexp, iexp), mks{iexp}, 'Color', cc(iexp,:),'MarkerSize', 10, 'MarkerFaceColor', cc(iexp,:), 'HandleVisibility', 'off')

end

plot(tmod, log10(LDL)*ones(nt, 1), '--', 'Color', 'k', 'LineWidth', 1.5)
text(11, log10(LDL) - 0.5, 'LDL', 'Interpreter', 'Latex', 'FontSize', 17)

xlim([0 48])
ylim([0 20])


set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'Latex', 'XTick', [0 12 24 36 48])

xlabel({'','Time (h)'}, 'Interpreter', 'Latex', 'FontSize', 17)

hold off
box off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Third subplot: fit results with process noise
clear logCFUS_data logCFUS_mod logCFUS_traj par

path_name = 'C:\Users\RYZEN9\Desktop\NEREA\Matlab\2023-HeteroRes\Paper_IdentObs\PE\TP_data';
file_name = strcat(path_name, filesep, 'resFit_TP_', num2str(Ntraj),'traj.mat');

load(file_name, 'optSol', 'nom_par', 'CFUS_data', 'CFUS_traj')

CFUS_traj(CFUS_traj < LDL) = LDL;
logCFUS_data = log10(CFUS_data);
logCFUS_traj = log10(CFUS_traj);

% Obtain optimal parameter values to plot model average:
optpar  = optSol.*nom_par;

bS        = optpar(1);
bR        = optpar(2);
alpha_b   = optpar(3);
d_maxS    = optpar(4);
beta_d    = optpar(5);
alpha_d   = optpar(6);
EC_50d    = optpar(7);
H_d       = optpar(8);
xi_SR     = optpar(9);
k_xi      = optpar(10);
N_T0      = optpar(11);
lambda_T0 = optpar(12);

% Initial condition:
f0  = exp(-lambda_T0*rr);
f0  = f0/sum(f0);
N0  = N_T0*f0;

% Initialice total population size:
logCFUS_mod = zeros(nt, Nexp);

% Mutation rates matrix:
Xi = xi_SR*exp(k_xi*(1 - RR));
Xi = Xi - diag(diag(Xi));

AA_aux = Xi' - diag(sum(Xi, 2));

% Calculate growth rate:
b = bS*bR./(bR + rr.^alpha_b*(bS - bR));

% Calculate death rate:
d_max = d_maxS*beta_d^alpha_d*(1 - rr.^alpha_d)./(beta_d^alpha_d + rr.^alpha_d);

% ----------------------------------------------------------------------- %
% Plot fit results with process noise:

ind_tplot = 1:nt;

transp = [0.6 0.6 0.6];

subplot(1, 3, 3)

hold on

for iexp = 1:Nexp
    
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
    [~, xout] = ode15s(@(t,s) Odes_cte(t, s, AA), tmod, N0, ODEoptions);
   
    % Total population:
    logCFUS_mod(1:nt, iexp) = log10(sum(xout, 2));
    
    for itraj = 1:Ntraj
        plot(tmod(ind_tplot), logCFUS_traj(ind_tplot, iexp, itraj), 'Color', [cc(iexp,:) transp(itraj)], 'HandleVisibility', 'off')
    end
    
    aux = reshape(CFUS_traj(1:nt, iexp, 1), nt, 1);
    
    NaN_ind = find(isnan(aux));
    
    logCFUS_mod(NaN_ind, iexp) = NaN;
    
    plot(tmod, logCFUS_mod(1:nt, iexp), 'Color', cc(iexp, :),'LineWidth', 2)
    
    plot(texp, logCFUS_data(1:ntexp, iexp), mks{iexp}, 'Color', cc(iexp,:),'MarkerSize', 10, 'MarkerFaceColor', cc(iexp,:), 'HandleVisibility', 'off')

end
plot(tmod, log10(LDL)*ones(nt, 1), '--', 'Color', 'k', 'LineWidth', 1.5)
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

indexOfInterest = find(tmod < xfz & tmod > x0z); % range of t near perturbation

aux_ind  = intersect(indexOfInterest, texp_ind);
aux_ind1 = find(texp_ind > indexOfInterest(1), 1 );
aux_ind2 = find(texp_ind < indexOfInterest(end), 1, 'last');


for iexp = 1:Nexp
    for itraj = 1:Ntraj
        plot(tmod(indexOfInterest),logCFUS_traj(indexOfInterest, iexp, itraj), 'Color', [cc(iexp,:) transp(itraj)], 'HandleVisibility', 'off') % plot on new axes
        hold on
        plot(tmod(indexOfInterest),logCFUS_mod(indexOfInterest, iexp),'Color', cc(iexp,:),'LineWidth', 2, 'HandleVisibility', 'off') % plot on new axes
        
    end
    if numel(aux_ind1) > 0
        plot(tmod(aux_ind), logCFUS_data(aux_ind1:aux_ind2, iexp), mks{iexp}, 'Color', cc(iexp,:),'MarkerSize', 10, 'MarkerFaceColor', cc(iexp,:), 'HandleVisibility', 'off')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save figure:
%set(gcf, 'Position', [-1473         281        1167         614])
set(gcf, 'Position', [136         137        1590         823])
print('-r720','PanelFitResHD', '-dpng')















