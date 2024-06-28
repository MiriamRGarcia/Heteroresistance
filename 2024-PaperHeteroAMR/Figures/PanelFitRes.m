%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PanelFitRes: Plot results of model calibration using MLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
%close all

% Add folder with necessary functions:
addpath('../Functions')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of replicates used for calibration:
m_traj    = 3;  

% Define number of run to load calibration results:
m_run     = 1;

% Number of AMR levels used for calibration (to load results):
m_r       = 50;

% Number of AMR levels used for generating data:
m_rdata   = 50;

% Minimum and maximum AMR levels used to generate data:
ra_data   = 0;
rb_data   = 1;

% Method used to generate BD trajectories for PN case ('SSA' or 'RSSA'):
method    = 'SSA';                                                           

% Low detection limit (only for plot considerations):
LDL       = 10; %(CFUS/mL)

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

%%
% Load common experimental setup:
file_name = sprintf('../Results/resPE/resPE_%usubpop_MNHo_%utraj_run%u.mat', m_r, m_traj, m_run);

load(file_name, 'r', 'tsim', 'texp', 'Cexp', 'ODEoptions')

% Indexes of experimental times within sampling times:
texp_ind = find(ismember(tsim, texp));

% Problem sizes:
m_t    = numel(tsim);
m_e    = numel(Cexp);
m_texp = numel(texp_ind);

r_data = linspace(ra_data, rb_data, m_rdata).';
R_data = repmat(r_data, 1, m_rdata) - repmat(r_data.', m_rdata, 1);                       
R_data = R_data - triu(R_data) + tril(R_data).';

R = repmat(r, 1, m_r) - repmat(r.', m_r, 1);                       
R = R - triu(R) + tril(R).';

% ----------------------------------------------------------------------- %
% First subplot (MNHo case):

% Generate random trajectories with exact parameter values:
% (note that this is saved in file in the sampling times texp, not in tsim,
% so it must be recalculated):
load(file_name, 'pars', 'seed', 'sd')

% Call to function simulating the average BD model:
if m_rdata > 2                                                     
    [~, N_Tdata] = Sim_aveBD(r_data, R_data, tsim, Cexp, pars, ODEoptions);
else
    [~, N_Tdata] = Sim_aveBD_2subpop(r_data, R_data, tsim, Cexp, pars, ODEoptions);
end


% Add measuremennt noise to average total counts:
rng(seed)                                                                  % Set seed for generating randoms;
               
N_Tdata     = log10(repmat(N_Tdata, 1, 1, m_traj)) + ...                   % Replicates used for calibration;
                      sd*randn(m_t, m_e, m_traj);
                  
N_Tave_data = sum(N_Tdata, 3)/m_traj;                                      % Average of the replicates used for calibration;
      

N_Tdata(N_Tdata < log10(LDL)) = log10(LDL);


% Load fit results:
load(file_name, 'pars_opt')

% Call function simulating average BD model with the optimal parameters:
if m_r > 2                                                     
    [~, N_Tmod] = Sim_aveBD(r, R, tsim, Cexp, pars_opt, ODEoptions);
else
    [~, N_Tmod] = Sim_aveBD_2subpop(r, R, tsim, Cexp, pars_opt, ODEoptions);
end

% Initialise legend:
lgd = cell(m_e, 1);

% Figure:
fig = figure;

set(gcf,'color','w');

transp    = [0.3 0.5 0.7];
ind_tplot = 1:1000:m_t;


% First subplot:
subplot(1, 3, 1)

hold on

for iexp = 1:m_e
    
    for itraj = 1:m_traj
        plot(tsim(ind_tplot), N_Tdata(ind_tplot, iexp, itraj), 'Color', [cc(iexp,:) transp(itraj)], 'HandleVisibility', 'off')
    end
    
    plot(tsim, log10(N_Tmod(1:m_t, iexp)), 'Color', cc(iexp, :),'LineWidth', 2)
    
    plot(texp, N_Tave_data(texp_ind, iexp), mks{iexp}, 'Color', cc(iexp,:),'MarkerSize', 10, 'MarkerFaceColor', cc(iexp,:), 'HandleVisibility', 'off')
    
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

%%
% ----------------------------------------------------------------------- %
% Second subplot (MNHe case):
file_name = sprintf('../Results/ResPE/resPE_%usubpop_MNHe_%utraj_run%u.mat', m_r, m_traj, m_run);

load(file_name, 'seed', 'var_a', 'var_b')

% Call to function simulating average BD model:
if m_rdata > 2                                                     
    [~, N_Tdata] = Sim_aveBD(r_data, R_data, tsim, Cexp, pars, ODEoptions);
else
    [~, N_Tdata] = Sim_aveBD_2subpop(r_data, R_data, tsim, Cexp, pars, ODEoptions);
end

% Generate random trajectories:
rng(seed)

var         = var_a*N_Tdata.^var_b;
N_Tdata     = repmat(N_Tdata, 1, 1, m_traj) + ...
                repmat(sqrt(var), 1, 1, m_traj).*randn(m_t, m_e, m_traj);
            
N_Tave_data = sum(N_Tdata, 3)/m_traj;

N_Tdata(N_Tdata < LDL)         = LDL;

% Load fit results:
load(file_name, 'pars_opt')

% Call function simulating average BD model with the optimal parameters:
if m_r > 2                                                     
    [~, N_Tmod] = Sim_aveBD(r, R, tsim, Cexp, pars_opt, ODEoptions);
else
    [~, N_Tmod] = Sim_aveBD_2subpop(r, R, tsim, Cexp, pars_opt, ODEoptions);
end

ind_tplot = 1:1000:m_t;

subplot(1, 3, 2)

hold on
    
for iexp = 1:m_e

    for itraj = 1:m_traj
        plot(tsim(ind_tplot), log10(N_Tdata(ind_tplot, iexp, itraj)), 'Color', [cc(iexp,:) transp(itraj)], 'HandleVisibility', 'off')
    end
    
     plot(tsim, log10(N_Tmod(1:m_t,iexp)), 'Color', cc(iexp, :),'LineWidth', 2)
     plot(texp, log10(N_Tave_data(texp_ind, iexp)), mks{iexp}, 'Color', cc(iexp,:),'MarkerSize', 10, 'MarkerFaceColor', cc(iexp,:), 'HandleVisibility', 'off')

end

plot(tsim, log10(LDL)*ones(m_t, 1), '--', 'Color', 'k', 'LineWidth', 1.5)
text(11, log10(LDL) - 0.5, 'LDL', 'Interpreter', 'Latex', 'FontSize', 17)

xlim([0 48])
ylim([0 20])


set(gca, 'FontSize', 15, 'TickLabelInterpreter', 'Latex', 'XTick', [0 12 24 36 48])

xlabel({'','Time (h)'}, 'Interpreter', 'Latex', 'FontSize', 17)

hold off
box off

%%
% ----------------------------------------------------------------------- %
% Third subplot (PN case):
file_name = sprintf('../Results/resPE/resPE_%usubpop_PN_%utraj_run%u.mat', m_r, m_traj, m_run);

load(file_name, 'pars_opt', 'itraj')

% Load SSA trajectories:
count = 1;

% Indices of the experiments to plot:
ind_Cexp = 1:5;
for ii = itraj
    file_name = sprintf('../Results/ResSSA/res%s_%03u', method, ii);
    load(file_name, 'N_T')
    
    N_Tdata(1:m_t, 1:m_e, count) = N_T(1:m_t, ind_Cexp);    
    count = count + 1;
end

N_Tave_data = sum(N_Tdata, 3)/m_traj;
N_Tave_data = N_Tave_data(texp_ind, 1:m_e);

N_Tdata(N_Tdata < LDL)         = LDL;
N_Tave_data(N_Tave_data < LDL) = LDL;

% Call function simulating average BD model with the optimal parameters:
if m_r > 2                                                     
    [~, N_Tmod] = Sim_aveBD(r, R, tsim, Cexp, pars_opt, ODEoptions);
else
    [~, N_Tmod] = Sim_aveBD_2subpop(r, R, tsim, Cexp, pars_opt, ODEoptions);
end

% ----------------------------------------------------------------------- %
% Plot fit results with process noise:

ind_tplot = 1:m_t;

transp = [0.6 0.6 0.6];

subplot(1, 3, 3)

hold on

for iexp = 1:m_e
    
    for itraj = 1:m_traj
        plot(tsim(ind_tplot), log10(N_Tdata(ind_tplot, iexp, itraj)), 'Color', [cc(iexp,:) transp(itraj)], 'HandleVisibility', 'off')
    end
    
    aux = reshape(N_Tdata(1:m_t, iexp, 1), m_t, 1);
    
    NaN_ind = find(isnan(aux));
    
    N_Tmod(NaN_ind, iexp) = NaN;
    
    plot(tsim, log10(N_Tmod(1:m_t, iexp)), 'Color', cc(iexp, :),'LineWidth', 2)
    
    plot(texp, log10(N_Tave_data(1:m_texp, iexp)), mks{iexp}, 'Color', cc(iexp,:),'MarkerSize', 10, 'MarkerFaceColor', cc(iexp,:), 'HandleVisibility', 'off')

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
        plot(tsim(indexOfInterest),log10(N_Tdata(indexOfInterest, iexp, itraj)), 'Color', [cc(iexp,:) transp(itraj)], 'HandleVisibility', 'off') % plot on new axes
        hold on
        plot(tsim(indexOfInterest),log10(N_Tmod(indexOfInterest, iexp)),'Color', cc(iexp,:),'LineWidth', 2, 'HandleVisibility', 'off') % plot on new axes
        
    end
    if numel(aux_ind1) > 0
        plot(tsim(aux_ind), log10(N_Tave_data(aux_ind1:aux_ind2, iexp)), mks{iexp}, 'Color', cc(iexp,:),'MarkerSize', 10, 'MarkerFaceColor', cc(iexp,:), 'HandleVisibility', 'off')
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

%%
% ----------------------------------------------------------------------- %
% Save figure if desired:
set(gcf, 'Position', [136         137        1590         823])

if fig_print == 1
    print('-r720', sprintf('PanelFitRes_%usubpop_%utraj_run%u', m_r, m_traj, m_run), fig_form)
end

rmpath('../Functions')
