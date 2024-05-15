%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MainPanelSens: Main file to generate panel of parameter sensitivities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

addpath('Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set precision to solve ODEs:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

% Save figure (=1) or not (=0):
fig_print = 0;

% Format to save figure:
fig_form  = '-dpng'; %'djpg';

% Parameter values:
b_S       = 0.63;
b_R       = 0.36;
alpha_b   = 2;

d_maxS    = 3.78;
beta_d    = 0.4;
alpha_d   = 3;

EC_50d    = 1;
H_d       = 1;

xi_SR     = 1e-6;
k_xi      = log(1e2);

N_T0      = 1e6;
lambda_T0 = 50;

pars      = [b_S;b_R;alpha_b;d_maxS;alpha_d;beta_d;EC_50d;H_d;xi_SR;k_xi;N_T0;lambda_T0];
    
np        = numel(pars);

% Time discretisation:
t0        = 0; 
tf        = 48;       
ht        = 1e-2;             
tmod      = t0:ht:tf;        
nt        = numel(tmod);    

% Discretisation of the AMR level:
ra        = 0;
rb        = 1;
nr        = 50;
r         = linspace(ra, rb, nr).';

% Antimicrobial concentrations (constant for each experiment):
MIC_S     = EC_50d*(b_S/(d_maxS - b_S))^(1/H_d);
C         = MIC_S*linspace(0, 8, 100);
nC        = numel(C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user-defined settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial population size:
f0 = exp(-lambda_T0*r);
f0 = f0/sum(f0);
N0 = floor(N_T0*f0);

% ----------------------------------------------------------------------- %
% Call to function calculating sensitivities:
[NT, sens_NT] = SensMultiExp(tmod, r, C, pars, N0, ODEoptions);

% ----------------------------------------------------------------------- %
% Plot sensitivities:

% Define parameter names:
name_pars = {'$b_S$','$b_R$','$\alpha_b$','$d_{maxS}$','$\alpha_d$',...
             '$\beta_d$','$EC_{50d}$','$H_d$','$\xi_{SR}$','$k_{\xi}$','$N_{T0}$','$\lambda_{T0}$'};


fig = figure;

set(gcf,'color','w');

fig.Position = [554 152 1128 783];

Cticks = linspace(C(1), C(nC), 4);

for ip = 1:np
    
    subplot(4, 3, ip)
    
    % ------------------------------------------ %
    % Calculate normalised sensitivities:
    sens_aux    = reshape(sens_NT(1:nt, ip, 1:nC), nt, nC);
    sc_sens_aux = pars(ip)*sens_aux./NT;
    
    % ------------------------------------------ %
    % Absolute value:
    sc_sens_aux = abs(sc_sens_aux);
    
    % ------------------------------------------ %
    % Contour plot:
    contourf(tmod, C, sc_sens_aux.', 100, 'edgecolor', 'none')
    
    % ------------------------------------------ %
    % Set scale for parameter N_T0:
    if ip == 11
        caxis([0 max(max(sc_sens_aux))]);
    end
    
    % ------------------------------------------ %
    % Colorbar of subplot:
    CB   = colorbar;
    lCB  = get(CB,'Limits');
    tCB  = linspace(lCB(1), lCB(2), 4);
    
    CB.FontSize = 11;
    
    set(CB,'Ticks',tCB,'TickLabelInterpreter', 'Latex')
    
    TLCB = arrayfun(@(x) sprintf('%.2f',x), tCB, 'un', 0);
    
    set(CB,'TickLabels',TLCB)
    
    % ------------------------------------------ %
    % General settings for subplot:    
    xlim([tmod(1) tmod(end)])
    ylim([C(1) C(end)])
    
    % Axis properties:
    set(gca,'XTick',[0 10 20 30 40 48], 'YTick', Cticks, 'TickLabelInterpreter', 'Latex', 'FontSize', 12)
    
    % Title:
    title(name_pars{ip}, 'Interpreter', 'Latex', 'FontSize', 16)
    
    % Set ticks on y axis:
    ax = gca;
    ax.YAxis.TickLabelFormat = '%.2f';
    
    box off
end

% ------------------------------------------ %
% Only one xy axes for the figure:
han = axes(fig,'visible','off'); 

han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';

Xlim = xlim;           
Xlb  = mean(Xlim); 
xlabel(han,{'','Time (h)'}, 'Interpreter','Latex','FontSize', 16, 'Position',[Xlb-0.03 -0.07],'HorizontalAlignment','center')
ylabel(han,{'$C$ (mg/L)',''}, 'Interpreter','Latex','FontSize', 16);

% ----------------------------------------------------------------------- %
% Print figure if desired:
if fig_print == 1
    print('-r720','PanelSens', fig_form)
end


rmpath('Functions')