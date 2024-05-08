%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot panel with simulation results of the time-kill assay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

addpath('Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) Load data of Gillespie trajectories:

% Results path:
path_name = 'C:\Users\RYZEN9\Desktop\NEREA\Matlab\2023-HeteroRes\2024-PaperHeteroAMR\SSA\Results';

% Load experimental setup:
file_name = 'resGill_001.mat';
file_name = strcat(path_name, filesep, file_name);

load(file_name, 'rr', 'tmod', 'Cexp', 'par')

% Problem sizes:
nr    = numel(rr);
nt    = numel(tmod);
Nexp  = numel(Cexp);
Ntraj = 1000;

% Times to almacenate the heteroresistance distribution:
AMR_tind   = 1:100:nt;
AMR_tt     = tmod(AMR_tind);
AMR_nt     = numel(AMR_tind);

% Initialise Gillespie trajectories of heteroresistance distribution (CFUS/mL):
CFUS_data  = zeros(AMR_nt, nr, Nexp, Ntraj);

% Initialise Gillespie trajectories of total cell counts (CFUS/mL):
CFUST_data = zeros(nt, Nexp, Ntraj);

% Initialise PDF of the AMR level:
pdf_AMR    = zeros(nt, nr, Nexp);

% Loop in the trajectories:
for itraj = 1:Ntraj

    % Create file name:
    file_name = sprintf('resGill_%03u', itraj);
    file_name = strcat(path_name, filesep, file_name);
    
    % Load results file:
    load(file_name, 'CFUS_exp', 'CFUST_exp')

    % Preprocesing data:
    for iexp = 1:Nexp
        aux       = [0;
                     CFUST_exp(1:(nt - 1), iexp)];
                 
        [rind, ~] = find((CFUST_exp(1:nt, iexp) - aux == 0) & (CFUST_exp(1:nt, iexp) - repmat(CFUST_exp(end, iexp), nt, 1) == 0) & (CFUST_exp(1:nt, iexp) - repmat(max(CFUST_exp(1:nt, iexp)), nt, 1) == 0));
        
        CFUST_exp(rind, iexp)      = NaN;
        
        CFUS_exp(rind, 1:nr, iexp) = NaN;
    end

    % Almacenate cell count data:
    CFUS_data(1:AMR_nt, 1:nr, 1:Nexp, itraj) = CFUS_exp(AMR_tind, 1:nr, 1:Nexp);
    
    % Almacenate total cell count data:  
    CFUST_data(1:nt, 1:Nexp, itraj) = CFUST_exp;
       
    % Calculate pdf of the AMR level:
    aux_CFUS     = CFUS_exp;
    aux_nT_CFUS  = repmat(CFUST_exp, nr, 1);
    
    aux_nT_CFUS  = reshape(aux_nT_CFUS, nt, nr, Nexp);

    pdf_AMR(1:nt, 1:nr, 1:Nexp) = pdf_AMR(1:nt, 1:nr, 1:Nexp) + aux_CFUS./aux_nT_CFUS;
 
end

% PDF of the AMR level:
pdf_AMR = pdf_AMR/Ntraj;

% Sample average of cell counts:
CFUS_ave_data  = sum(CFUS_data, 4)/Ntraj;

% Sample average of total cell counts:
CFUST_ave_data = sum(CFUST_data, 3)/Ntraj;

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
f0  = exp(-lambda_T0*rr);
f0  = f0/sum(f0);
N0  = N_T0*f0;

% Initialice population sizes:
CFUS_mod  = zeros(nt, nr, Nexp);
CFUST_mod = zeros(nt, Nexp);

% Initialise coefficient matrix:
RR = repmat(rr, 1, nr) - repmat(rr.', nr, 1);                      
RR = RR - triu(RR) + tril(RR).';

Xi = xi_SR*exp(k_xi*(1 - RR));
Xi = Xi - diag(diag(Xi));

AA_aux = Xi.' - diag(sum(Xi, 2));

% Calculate growth rate at current time:
b = bS*bR./(bR + rr.^alpha_b*(bS - bR));

% Calculate kill rate at current time:
d_max = d_maxS*beta_d^alpha_d*(1 - rr.^alpha_d)./(beta_d^alpha_d + rr.^alpha_d);

for iexp = 1:Nexp
    
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
    [~, xout] = ode15s(@(t, s) Odes_cte(t, s, AA), tmod, N0, ODEoptions);
    
    %-------------------------------------------------%
    % Almacenate cell numbers for current antimicrobial concentration:
    CFUS_mod(1:nt, 1:nr, iexp) = xout;
    CFUST_mod(1:nt, iexp)      = sum(xout, 2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (4) Plot results (fitness cost and kill rate):
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

set(fig, 'Color', 'w')
set(gcf, 'Position', [173 99 1586 882])

nr_plot = 1e3;
rr_plot = linspace(0, 1, nr_plot);

% Calculate growth rate:
b_plot = bS*bR./(bR + rr_plot.^alpha_b*(bS - bR));

% ----------------------------------------------------------------------- %
% Plot birth-rate:
ss = subplot(2, 3, 1);
hold on
plot(rr_plot, b_plot, 'Color', bb5, 'LineWidth', 1.5)
plot(rr, bS*ones(nr, 1), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
plot(rr, bR*ones(nr, 1), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
text(0.1, bS + 0.02, '$b_S$', 'Interpreter', 'Latex', 'FontSize', 17)
text(0.1, bR + 0.02, '$b_R$', 'Interpreter', 'Latex', 'FontSize', 17)

rr_max = (bR*(alpha_b - 1)/((alpha_b + 1)*(bS - bR)))^(1/alpha_b);
plot(rr_max*ones(nr, 1), linspace(bR, bS, nr), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
ttb = text(rr_max, 0.32, '$r_b$', 'Interpreter', 'Latex', 'FontSize', 17);
ttb.Position = [0.661498708010336 0.333680781758958 0];

hold off

xlabel('$r$', 'Interpreter', 'Latex', 'FontSize', 17)
ylabel('Birth rate $b(r)$ (1/h)', 'Interpreter', 'Latex', 'FontSize', 17)

xticks = [0 0.2 0.4 0.6 0.8 1.0];
set(ss, 'XTick', xticks, 'TickLabelInterpreter', 'Latex', 'FontSize', 15)
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
plot(rr, d_maxS*ones(nr, 1), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
text(0.1, d_maxS + 0.3, '$d_{maxS}$', 'Interpreter', 'Latex', 'FontSize', 17)

rr_max = ((alpha_d - 1)*beta_d^alpha_d/(alpha_d + 1))^(1/alpha_d);
plot(rr_max*ones(nr, 1), linspace(0, d_maxS, nr), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5)
ttk = text(rr_max, -0.1, '$r_d$', 'Interpreter', 'Latex', 'FontSize', 17);
ttk.Position = [0.286472458455655 -0.217263843648209 0];
hold off

xlabel('$r$', 'Interpreter', 'Latex', 'FontSize', 17)
ylabel('Death rate $d(r,C)$ (1/h)', 'Interpreter', 'Latex', 'FontSize', 17)

set(ss, 'XTick', xticks, 'TickLabelInterpreter', 'Latex','FontSize', 15)
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

for iexp = Exps
    
    NaNind = find(isnan(logCFUST_ave_data(1:nt, iexp)));
    
    % ----------------------------------------------- %
    % Plot model average of data:
    logCFUST_mod(NaNind, iexp) = NaN;
    
    plot(tmod, logCFUST_mod(1:nt, iexp), 'Color', cc(iexp, :), 'LineWidth', 1.0)
    
    % ----------------------------------------------- %
    % Plot average value of trajectories:
    plot(tmod, logCFUST_ave_data(1:nt, iexp), 'LineStyle', '--', 'Color', cc(iexp, :), 'LineWidth', 1.0, 'HandleVisibility', 'off')

    % ----------------------------------------------- %
    % Shade area between quartile 0.25 and 0.75:
    lowlim = logCFUST_data_quarlow(1:nt, iexp);
    uplim  = logCFUST_data_quarup(1:nt, iexp);
    
    if numel(NaNind) > 0
        lowlim(NaNind:end) = [];
        uplim(NaNind:end)  = [];
    
        taux = tmod(1:NaNind - 1);
    else
        taux = tmod;
    end
    
    coord_up  = [taux.', uplim];
    coord_low = [taux.', lowlim];
    
    coord_combine = [coord_up;flipud(coord_low)];
    
    fill(coord_combine(:,1),coord_combine(:,2), cc(iexp, :), 'EdgeColor', 'None', 'FaceAlpha', 0.4, 'HandleVisibility', 'off')
    

    % ----------------------------------------------- %
    % Maximum an minimum of data without outliers:
    CFUST_data_aux = reshape(CFUST_data(1:nt, iexp, 1:Ntraj), nt, Ntraj);
    IQR_aux        = IQR(1:nt, iexp);
    limlow         = repmat(CFUST_data_quarlow(1:nt, iexp) - 1.5*IQR_aux, 1, Ntraj);
    limup          = repmat(CFUST_data_quarup(1:nt, iexp)  + 1.5*IQR_aux, 1, Ntraj);
    outind         = find(limup < CFUST_data_aux | CFUST_data_aux < limlow);
     
    % ----------------------------------------------- %
    % Scatter plot with data outliers:
    %taux2 = repmat(tmod.', 1, Ntraj);
    
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
    lgd{lgd_count} = sprintf('$C=%.1f$ mg/L', Cexp(iexp));%strcat('$C=\;$ ', num2str(Cexp(iexp)), ' (mg/L)');
    
    lgd_count      = lgd_count + 1;

end

% Low detection limit:
LDL = log10(10);

plot(tmod, LDL*ones(nt,1), 'k', 'LineStyle', '--', 'LineWidth', 1.5)
xlabel('Time (h)', 'Interpreter', 'Latex', 'FontSize', 17)
ylabel('Total cell count ($log_{10}$CFUS/mL)', 'Interpreter', 'Latex', 'FontSize', 17)

xlim([0 40])

set(gca, 'XTick', [0.0 5.0 10.0 15.0 20.0 25.0 30.0 35.0 40.0])

legend(lgd, 'Orientation', 'Horizontal', 'Interpreter', 'Latex', 'edgecolor', 'none', 'FontSize', 17, 'Position',[0.183997965068409 0.00912237826224061 0.676653187471205 0.0297365115607248])
text(10, LDL + 0.1, 'LDL', 'Interpreter', 'Latex', 'FontSize', 17)

hold off
box off
set(ss, 'TickLabelInterpreter', 'Latex', 'FontSize', 15)
%set(ss, 'FontSize', 13)
%ss.XAxis.TickLabelFormat = '%,.1f';
%ss.YAxis.TickLabelFormat = '%,.1f';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (6) Plot results Waterfall (pdf of the AMR level at different times):

% Index of intersection between experiments:
ind45i = find(CFUST_ave_data(:,5) > CFUST_ave_data(:,4), 1);
ind45f = find(CFUST_ave_data(:,5) > CFUST_ave_data(:,4), 1, 'last');
ind49i = find(CFUST_ave_data(:,9) > CFUST_ave_data(:,4), 1);
ind49f = find(CFUST_ave_data(:,9) > CFUST_ave_data(:,4), 1, 'last');
ind59i = find(CFUST_ave_data(:,9) > CFUST_ave_data(:,5), 1);
ind59f = find(CFUST_ave_data(:,9) > CFUST_ave_data(:,5), 1, 'last');

tt1 = find(tmod >= 5.6, 1);
tt2 = find(tmod >= 10, 1);
tt3 = find(tmod >= 30, 1);
%tt = find(tmod ==[5 9 13 17 21 25 48]);

tind_pdfAMR = sort([ind49i ind49f ind45i ind45f ind59i ind59f 30001]);
%tind_pdfAMR = [ind59i ind49i 14001 20001 24001 30001];

nt_pdfAMR   = numel(tind_pdfAMR);

% Experiments to plot:
Exps        = [4 5 9];
Nexp_aux    = numel(Exps);

% Set baseline and basestep:
basestep    = 0.15;
baseline    = zeros(nr, 1) + basestep*(nt_pdfAMR - 1);
%baseline    = zeros(nr, 1);

% Abcises to plot text with times:
rtext  = rr([6 8 10 12 14 16 18]);%rr(10)*ones(1,nt_pdfAMR);

% Transparency for the experiments:
transp = [1.0 0.6 0.6];


ss = subplot(2, 3, [3 6]);

hold on

% For loop on times to plot:
for it = 1:nt_pdfAMR
    
    plot_count = 1;
    it         = nt_pdfAMR - it + 1;
    indt       = tind_pdfAMR(it);
    
    % Set axis to calculate area above:
    lowlim     = baseline;
    
    baseline   = baseline - basestep;%baseline + basestep;

    for iexp = Exps
        
        NaNind = find(isnan(pdf_AMR(indt, 1, iexp)));
        
        if numel(NaNind) < 1
            
            % ----------------------------------------------- %
            % Plot distribution:      
            uplim         = pdf_AMR(indt, 1:nr, iexp).' + lowlim;
    
            coord_up      = [rr, uplim];
            coord_low     = [rr, lowlim];
    
            coord_combine = [coord_up;flipud(coord_low)];
            
            fill(coord_combine(:,1),coord_combine(:,2), cc(iexp, :), 'EdgeColor', bb5, 'FaceAlpha', transp(plot_count))
        end
        
        plot_count = plot_count + 1;
    end
    
    % ----------------------------------------------- %
    % Put text label with time:    
    ytext = lowlim(1) + 0.03;
    text(rtext(it), ytext, sprintf('$t=%.1f$h', tmod(indt)), 'Interpreter', 'Latex','FontSize', 17)
     
end

yticks = 0:basestep:basestep*(nt_pdfAMR);
set(ss, 'XTick', xticks, 'YTick', yticks, 'TickLabelInterpreter', 'Latex','FontSize',15)
%set(ss, 'YTick', yticks, 'FontSize',13)
%set(ss, 'TickLabelInterpreter', 'Latex')
%ss.XAxis.TickLabelFormat = '%,.1f';
ss.YAxis.TickLabelFormat = '%,.2f';

% Add grid:
grid on
ss.YMinorGrid         = 'on';
ss.MinorGridColor     = cc(5,:);
ss.MinorGridLineStyle = '--';
ss.MinorGridAlpha     = 0.3;
ss.GridColor          = cc(5,:);
ss.GridLineStyle      = '--';
ss.GridAlpha          = 0.3;
ss.Layer              = 'top';

xlabel('$r$', 'Interpreter', 'Latex','FontSize', 17)
ylabel('Subpopulation frequency', 'Interpreter', 'Latex','FontSize', 17)
ylim([0 1.05])


print('-r720', 'PanelSimResHD', '-dpng')

rmpath('Functions')