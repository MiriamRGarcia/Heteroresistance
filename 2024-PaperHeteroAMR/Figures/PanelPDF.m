%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PanelPDF: Plot panel with heteroresistance distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

addpath('../Functions')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined options:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of trajectories:
m_traj = 1000;

% Choose implementation (direct method = SSA or rejection based = RSSA):
method = 'SSA'; % = 'SSA'; = 'RSSA';

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

%%
% ----------------------------------------------------------------------- %
% Load data of Gillespie trajectories:

% Results path:
file_name = sprintf('../Results/ResSSA/res%s_001.mat', method);
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
    file_name = sprintf('../Results/ResSSA/res%s_%03u', method, itraj);

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

mt_pdfAMR   = numel(tind_pdfAMR);


% ----------------------------------------------------------------------- %
% Options for area plot:
ra = r(1);
rb = r(m_r);

basestep = 0.15;                                                           % Set step between baselines;
baseline = zeros(m_r, 1) + basestep*(mt_pdfAMR - 1);                       % Set base lines for area plots;
rtext    = linspace(ra + 0.1, rb - 0.1, mt_pdfAMR);%r([6 8 10 12 14 16]);  % x position of text in figure;
transp   = [1.0 0.6 0.6];                                                  % Transparency for the area plots;
xticks   = [0 0.2 0.4 0.6 0.8 1.0];

fig = figure;

set(gcf, 'Position', [237.6667  138.3333  972.6667  479.6667], 'Color', 'w')

ss = subplot(1, 3, 1);

hold on

% For loop on times to plot:
for it = 1:mt_pdfAMR
    
    plot_count = 1;
    aux_it     = mt_pdfAMR - it + 1;
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

yticks = 0:basestep:basestep*(mt_pdfAMR);
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
R  = repmat(r, 1, m_r) - repmat(r.', m_r, 1);                      
R  = R - triu(R) + tril(R).';

if m_r > 2
    [N_mod, N_Tmod] = Sim_aveBD(r, R, tsim, Cexp, pars, ODEoptions);
else
    [N_mod, N_Tmod] = Sim_aveBD_SR(r, R, tsim, Cexp, pars, ODEoptions);
end

% ----------------------------------------------------------------------- %
% Plot results:

% Experiments to plot:
Exps        = [5 9];

tind_pdfAMR    = {};
tind_pdfAMR{1} = sort([ind49i ind49f ind45i ind45f ind59i ind59f]);
tind_pdfAMR{2} = sort([ind49i ind49f ind45i ind45f ind59i ind59f]);

mt_pdfAMR   = numel(tind_pdfAMR{1});


% Set baseline:
baseline    = zeros(m_r, 1);

ccaux(1,:)  = [169,188,173]/256;%[47,126,68]/256;
cccaux(1,:) = [30,69,52]/256;
ccaux(2,:)  = [216,203,175]/256;%[233,144,20]/256;%[175,114,0]/256;
cccaux(2,:) = [233,144,20]/256;

plot_count  = 1;

clear rtext
rtext{1}    = linspace(ra + 0.1, rb - 0.1, mt_pdfAMR);
rtext{2}    = linspace(ra + 0.1, rb - 0.1, mt_pdfAMR);

for it = 1:mt_pdfAMR
    
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
            uplim         = N_mod(indt, 1:m_r, iexp).'/N_Tmod(indt, iexp) + lowlim;
            
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
ss.GridColor          = cc(5,:);
ss.GridLineStyle      = '--';
ss.GridAlpha          = 0.3;
ss.Layer              = 'top';

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
ss.GridColor          = cc(9,:);
ss.GridLineStyle      = '--';
ss.GridAlpha          = 0.3;


if fig_print == 1
    print('-r720','PanelPDF', fig_form)
end

rmpath('../Functions')
