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
m_traj = 1000;

% Choose implementation (direct method = SSA or rejection based = RSSA):
method = 'RSSA'; % = 'SSA'; = 'RSSA';

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
file_name = fprintf('../SSA/Results/res%s_001.mat', method);

load(file_name, 'r', 'tsim', 'Cexp', 'pars')

% Problem sizes:
m_r = numel(r);
m_t = numel(tsim);
m_e = numel(Cexp);

% Initialise Gillespie trajectories of total cell counts (CFUS/mL):
N_Tdata = zeros(m_t, m_e, m_traj);

% Initialise PDF of the AMR level:
pdf_AMR = zeros(m_t, m_r, m_e);

% Loop in the trajectories:
for itraj = 1:m_traj

    % Create file name:
    file_name = sprintf('../SSA/Results/res%s_%03u', method, itraj);

    % Load results file:
    load(file_name, 'N', 'N_T')

    % Almacenate total cell count data:  
    N_Tdata(1:m_t, 1:m_e, itraj) = N_T;
       
    % Calculate pdf of the AMR level:
    aux_N   = N;
    aux_N_T = repmat(N_T, m_r, 1);
    
    aux_N_T = reshape(aux_N_T, m_t, m_r, m_e);

    pdf_AMR(1:m_t, 1:m_r, 1:m_e) = pdf_AMR(1:m_t, 1:m_r, 1:m_e) + aux_N./aux_N_T;
 
end

% PDF of the AMR level:
pdf_AMR     = pdf_AMR/m_traj;

% Sample average of total cell counts:
N_Tave_data = sum(N_Tdata, 3)/m_traj;


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
baseline = zeros(m_r, 1) + basestep*(nt_pdfAMR - 1);                        % Set base lines for area plots;
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
            uplim         = pdf_AMR(indt, 1:m_r, iexp).' + lowlim;
    
            coord_up      = [r, uplim];
            coord_low     = [r, lowlim];
    
            coord_combine = [coord_up;flipud(coord_low)];
            
            fill(coord_combine(:,1),coord_combine(:,2), cc(iexp, :), 'EdgeColor', bb, 'FaceAlpha', transp(plot_count))
        end
        
        plot_count = plot_count + 1;
        
    end
       
    ytext = lowlim(1) + 0.03;
    text(rtext(aux_it), ytext, sprintf('$t=%.1f$h', tsim(indt)), 'Interpreter', 'Latex','FontSize', 12)
     
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
Nmod  = zeros(m_t, m_r, m_e);
N_Tmod = zeros(m_t, m_e);

% Initialise coefficient matrix:
R  = repmat(r, 1, m_r) - repmat(r.', m_r, 1);                      
R  = R - triu(R) + tril(R).';

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
    Nmod(1:m_t, 1:m_r, iexp) = xout;
    N_Tmod(1:m_t, iexp)     = sum(xout, 2);

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
baseline    = zeros(m_r, 1);

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
            uplim         = pdf_AMR(indt, 1:m_r, iexp).' + lowlim;
    
            coord_up      = [r, uplim];
            coord_low     = [r, lowlim];
    
            coord_combine = [coord_up;flipud(coord_low)];
            
            subplot(1, 3, plot_count + 1)
            
            hold on
            
            fill(coord_combine(:,1),coord_combine(:,2), cc(iexp, :), 'EdgeColor', cccaux(plot_count,:), 'FaceAlpha', 1.0)
            
            % ----------------------------------------------- %
            % Plot model distribution:      
            uplim         = Nmod(indt, 1:m_r, iexp).'/N_Tmod(indt, iexp) + lowlim;
            
            coord_up      = [r, uplim];
    
            coord_combine = [coord_up;flipud(coord_low)];
            
            fill(coord_combine(:,1),coord_combine(:,2), ccaux(plot_count, :), 'EdgeColor', cccaux(plot_count,:), 'FaceAlpha', 0.5)
            
            % ----------------------------------------------- %
            % Put text label with time:    
            ytext = lowlim(1) + 0.03;
            text(rtext{plot_count}(it), ytext, sprintf('$t=%.1f$h', tsim(indt)), 'Interpreter', 'Latex','FontSize', 12)
            
            hold off
            
            plot_count = plot_count + 1;
        end
    end
    
    plot_count = 1;
end
ss = subplot(1, 3, 2);
yticks = 0:0.15:0.9;
set(ss, 'XTick', xticks, 'YTick', yticks, 'TickLabelInterpreter', 'Latex', 'FontSize', 11)
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

xlabel('$r$', 'Interpreter', 'Latex','FontSize', 12)


ss = subplot(1, 3, 3);
set(ss, 'XTick', xticks, 'YTick', yticks, 'TickLabelInterpreter', 'Latex', 'FontSize', 11)
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


if fig_print == 1
    print('-r720','PanelPDF', fig_form)
end
