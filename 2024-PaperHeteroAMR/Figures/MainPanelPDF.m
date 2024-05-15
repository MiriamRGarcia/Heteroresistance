%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MainPanelPDF: Plot panel with heteroresistance distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all
clc

addpath('Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined options:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of trajectories:
Ntraj = 1000;

% Define colors of plot:
bb   = [132,114,68]/256;

cc1  = [39, 183, 222]/256;
cc2  = [76,57,87]/256;
cc3  = [129,23,27]/256;                                            
cc4  = [237,106,90]/256;
cc5  = [51,115,87]/256;
cc6  = [250,163,0]/256;
cc7  = [250,163,0]/256;
cc9  = [250,163,0]/256;
cc8  = [250,163,0]/256;
cc10 = [250,163,0]/256;

cc   = [cc1;cc2;cc3;cc4;cc5;cc6;cc7;cc8;cc9;cc10];

% Save figure (=1) or not (=0):
fig_print = 0;

% Format to save figure:
fig_form  = '-dpng'; %'djpg';

% Set ODE solver precision:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user-defined settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------------------------------------------------------------------- %
% Load data of Gillespie trajectories:

% Results path:
file_name = '../SSA/Results/resSSA_001.mat';

load(file_name, 'r', 'tmod', 'Cexp', 'pars')

% Problem sizes:
nr   = numel(r);
nt   = numel(tmod);
Nexp = numel(Cexp);

% Initialise Gillespie trajectories of total cell counts (CFUS/mL):
N_Tdata = zeros(nt, Nexp, Ntraj);

% Initialise PDF of the AMR level:
pdf_AMR = zeros(nt, nr, Nexp);

% Loop in the trajectories:
for itraj = 1:Ntraj

    % Create file name:
    file_name = sprintf('../SSA/Results/resSSA_%03u.mat', itraj);

    % Load results file:
    load(file_name, 'N', 'N_T')
    
    % Preprocesing data:
    for iexp = 1:Nexp
        aux       = [0;
                     N_T(1:(nt - 1), iexp)];
                 
        [rind, ~] = find((N_T(1:nt, iexp) - aux == 0) & (N_T(1:nt, iexp) - repmat(N_T(end, iexp), nt, 1) == 0) & (N_T(1:nt, iexp) - repmat(max(N_T(1:nt, iexp)), nt, 1) == 0));
        
        N_T(rind, iexp)     = NaN;
        
        N(rind, 1:nr, iexp) = NaN;
    end

    % Almacenate total cell count data:  
    N_Tdata(1:nt, 1:Nexp, itraj) = N_T;
       
    % Calculate pdf of the AMR level:
    aux_N   = N;
    aux_N_T = repmat(N_T, nr, 1);
    
    aux_N_T = reshape(aux_N_T, nt, nr, Nexp);

    pdf_AMR(1:nt, 1:nr, 1:Nexp) = pdf_AMR(1:nt, 1:nr, 1:Nexp) + aux_N./aux_N_T;
 
end

% PDF of the AMR level:
pdf_AMR     = pdf_AMR/Ntraj;

% Sample average of total cell counts:
N_Tave_data = sum(N_Tdata, 3)/Ntraj;


% ----------------------------------------------------------------------- %
% Plot heteroresistance distribution at different times and experiments:

% Experiments to plot:
Exps     = [4 5 9];
Nexp_aux = numel(Exps);

% Times to plot (times of intersection between experiments):
ind45i = find(N_Tave_data(:,5) > N_Tave_data(:,4), 1);
ind45f = find(N_Tave_data(:,5) > N_Tave_data(:,4), 1, 'last');
ind49i = find(N_Tave_data(:,9) > N_Tave_data(:,4), 1);
ind49f = find(N_Tave_data(:,9) > N_Tave_data(:,4), 1, 'last');
ind59i = find(N_Tave_data(:,9) > N_Tave_data(:,5), 1);
ind59f = find(N_Tave_data(:,9) > N_Tave_data(:,5), 1, 'last');

tind_pdfAMR = sort([ind49i ind49f ind45i ind45f ind59i ind59f]);

nt_pdfAMR   = numel(tind_pdfAMR);


% ----------------------------------------------------------------------- %
% Options for area plot:

basestep = 0.15;                                                           % Set step between baselines;
baseline = zeros(nr, 1) + basestep*(nt_pdfAMR - 1);                        % Set base lines for area plots;
rtext    = r([6 8 10 12 14 16]);                                          % x position of text in figure;
transp   = [1.0 0.6 0.6];                                                  % Transparency for the area plots;
xticks   = [0 0.2 0.4 0.6 0.8 1.0];

fig = figure;

set(gcf, 'Position', [237.6667  138.3333  972.6667  479.6667], 'Color', 'w')

ss = subplot(1, 3, 1);

hold on

% For loop on times to plot:
for it = 1:nt_pdfAMR
    
    plot_count = 1;
    aux_it     = nt_pdfAMR - it + 1;
    indt       = tind_pdfAMR(aux_it);

    lowlim     = baseline;   
    baseline   = baseline - basestep;

    for iexp = Exps
        
        NaNind = find(isnan(pdf_AMR(indt, 1, iexp)));
        
        if numel(NaNind) < 1
            
            % ----------------------------------------------- %
            % Plot distribution:      
            uplim         = pdf_AMR(indt, 1:nr, iexp).' + lowlim;
    
            coord_up      = [r, uplim];
            coord_low     = [r, lowlim];
    
            coord_combine = [coord_up;flipud(coord_low)];
            
            fill(coord_combine(:,1),coord_combine(:,2), cc(iexp, :), 'EdgeColor', bb, 'FaceAlpha', transp(plot_count))
        end
        
        plot_count = plot_count + 1;
        
    end
       
    ytext = lowlim(1) + 0.03;
    text(rtext(aux_it), ytext, sprintf('$t=%.1f$h', tmod(indt)), 'Interpreter', 'Latex','FontSize', 12)
     
end

yticks = 0:basestep:basestep*(nt_pdfAMR);
set(ss, 'XTick', xticks, 'YTick', yticks, 'TickLabelInterpreter', 'Latex','FontSize',11)
ss.YAxis.TickLabelFormat = '%,.2f';

% Add minor grid:
grid on
ss.YMinorGrid         = 'on';
ss.MinorGridColor     = cc(5,:);
ss.MinorGridLineStyle = '--';
ss.MinorGridAlpha     = 0.3;
ss.GridColor          = cc(5,:);
ss.GridLineStyle      = '--';
ss.GridAlpha          = 0.3;
ss.Layer              = 'top';

xlabel('$r$', 'Interpreter', 'Latex','FontSize', 12)
%ylabel('Subpopulation frequency', 'Interpreter', 'Latex','FontSize', 12)
ylim([0 0.9])

% ----------------------------------------------------------------------- %
% Plot of deterministic VS stochastic heteroresistance distribution:

% Calculate deterministic average:

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
Nmod  = zeros(nt, nr, Nexp);
N_Tmod = zeros(nt, Nexp);

% Initialise coefficient matrix:
R  = repmat(r, 1, nr) - repmat(r.', nr, 1);                      
R  = R - triu(R) + tril(R).';

Xi = xi_SR*exp(k_xi*(1 - R));
Xi = Xi - diag(diag(Xi));

AA_aux = Xi.' - diag(sum(Xi, 2));

% Calculate growth rate at current time:
b = b_S*b_R./(b_R + r.^alpha_b*(b_S - b_R));

% Calculate kill rate at current time:
d_max = d_maxS*beta_d^alpha_d*(1 - r.^alpha_d)./(beta_d^alpha_d + r.^alpha_d);

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
    [~, xout] = ode15s(@(t, s) Odes_cte(t, s, AA), tmod, N_0, ODEoptions);
    
    %-------------------------------------------------%
    % Almacenate cell numbers for current antimicrobial concentration:
    Nmod(1:nt, 1:nr, iexp) = xout;
    N_Tmod(1:nt, iexp)     = sum(xout, 2);

end

% ----------------------------------------------------------------------- %
% Plot results:

% Experiments to plot:
Exps        = [5 9];

tind_pdfAMR    = {};
tind_pdfAMR{1} = sort([ind49i ind49f ind45i ind45f ind59i ind59f]);
tind_pdfAMR{2} = sort([ind49i ind49f ind45i ind45f ind59i ind59f]);

nt_pdfAMR   = numel(tind_pdfAMR{1});


% Set baseline:
baseline    = zeros(nr, 1);

ccaux(1,:)  = [169,188,173]/256;%[47,126,68]/256;
cccaux(1,:) = [30,69,52]/256;
ccaux(2,:)  = [216,203,175]/256;%[233,144,20]/256;%[175,114,0]/256;
cccaux(2,:) = [233,144,20]/256;

plot_count  = 1;

clear rtext
rtext{1}    = r([6 8 10 12 14 16]);
rtext{2}    = r([6 8 10 12 14 16]);

for it = 1:nt_pdfAMR
    
    lowlim   = baseline;

    baseline = baseline + basestep;  
    
    for iexp = Exps
        indt   = tind_pdfAMR{plot_count}(it);
        
        NaNind = find(isnan(pdf_AMR(indt, 1, iexp)));
        
        if numel(NaNind) < 1
            
            % ----------------------------------------------- %
            % Plot real distribution:      
            uplim         = pdf_AMR(indt, 1:nr, iexp).' + lowlim;
    
            coord_up      = [r, uplim];
            coord_low     = [r, lowlim];
    
            coord_combine = [coord_up;flipud(coord_low)];
            
            subplot(1, 3, plot_count + 1)
            
            hold on
            
            fill(coord_combine(:,1),coord_combine(:,2), cc(iexp, :), 'EdgeColor', cccaux(plot_count,:), 'FaceAlpha', 1.0)
            
            % ----------------------------------------------- %
            % Plot model distribution:      
            uplim         = Nmod(indt, 1:nr, iexp).'/N_Tmod(indt, iexp) + lowlim;
            
            coord_up      = [r, uplim];
    
            coord_combine = [coord_up;flipud(coord_low)];
            
            fill(coord_combine(:,1),coord_combine(:,2), ccaux(plot_count, :), 'EdgeColor', cccaux(plot_count,:), 'FaceAlpha', 0.5)
            
            % ----------------------------------------------- %
            % Put text label with time:    
            ytext = lowlim(1) + 0.03;
            text(rtext{plot_count}(it), ytext, sprintf('$t=%.1f$h', tmod(indt)), 'Interpreter', 'Latex','FontSize', 12)
            
            hold off
            
            %grid on
            
            plot_count = plot_count + 1;
        end
    end
    
    plot_count = 1;
end
ss = subplot(1, 3, 2);
yticks = 0:0.15:0.9;
set(ss, 'XTick', xticks, 'YTick', yticks, 'TickLabelInterpreter', 'Latex', 'FontSize', 11)
%ss.XAxis.TickLabelFormat = '%,.1f';
ss.YAxis.TickLabelFormat = '%,.2f';

% Add grid:
grid on
ss.YMinorGrid = 'on';
ss.MinorGridColor     = cc(5,:);
ss.MinorGridLineStyle = '--';
ss.MinorGridAlpha     = 0.3;
ss.GridColor     = cc(5,:);
ss.GridLineStyle = '--';
ss.GridAlpha     = 0.3;
ss.Layer         = 'top';
%box on

xlabel('$r$', 'Interpreter', 'Latex','FontSize', 12)


ss = subplot(1, 3, 3);
set(ss, 'XTick', xticks, 'YTick', yticks, 'TickLabelInterpreter', 'Latex', 'FontSize', 11)
%ss.XAxis.TickLabelFormat = '%,.1f';
ss.YAxis.TickLabelFormat = '%,.2f';
xlabel('$r$', 'Interpreter', 'Latex','FontSize', 12)

% Add grid:
grid on
ss.YMinorGrid = 'on';
ss.MinorGridColor     = cc(9,:);
ss.MinorGridLineStyle = '--';
ss.MinorGridAlpha     = 0.3;
ss.GridColor     = cc(9,:);
ss.GridLineStyle = '--';
ss.GridAlpha     = 0.3;
%ss.Layer         = 'top';
%box on

% Define ylabels:
% indt    = tind_pdfAMR{1};
% yticks  = [0 0.1 0.2 0.3 0.4 0.5];
% ylabels = {sprintf('$t=%.2f$', tmod(indt(1))),sprintf('$t=%.2f$', tmod(indt(2))),sprintf('$t=%.2f$', tmod(indt(3))),...
%     sprintf('$t=%.2f$', tmod(indt(4))),sprintf('$t=%.2f$', tmod(indt(5))),sprintf('$t=%.2f$', tmod(indt(6)))};
% 
% ss = subplot(1, 2, 1);
% set(ss, 'YTick', yticks, 'YTickLabel', ylabels, 'TickLabelInterpreter', 'Latex')
% ylim([0 0.7])
% 
% indt    = tind_pdfAMR{2};
% ylabels = {sprintf('$t=%.2f$', tmod(indt(1))),sprintf('$t=%.2f$', tmod(indt(2))),sprintf('$t=%.2f$', tmod(indt(3))),...
%     sprintf('$t=%.2f$', tmod(indt(4))),sprintf('$t=%.2f$', tmod(indt(5))),sprintf('$t=%.2f$', tmod(indt(6)))};
% 
% ss = subplot(1, 2, 2);
% set(ss, 'YTick', yticks, 'YTickLabel', ylabels, 'TickLabelInterpreter', 'Latex')
% ylim([0 0.7])


if fig_print == 1
    print('-r720','PanelPDF', fig_form)
end





