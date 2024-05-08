%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot panel of heteroresistance distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

addpath('Functions')

% Define colors of plot:
bb = [132,114,68]/256;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Load data of Gillespie trajectories:

% Number of trajectories:
Ntraj = 1000;

% Results path:
path_name = 'C:\Users\Usuario\Desktop\2024-PaperHeteroAMR\SSA\Results';

% Load experimental setup:
file_name = 'resGill_001.mat';
file_name = strcat(path_name, filesep, file_name);

load(file_name, 'rr', 'tmod', 'Cexp', 'par')

% Problem sizes:
nr    = numel(rr);
nt    = numel(tmod);
Nexp  = numel(Cexp);

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

% Sample average of total cell counts:
CFUST_ave_data = sum(CFUST_data, 3)/Ntraj;

%%
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Plot heteroresistance distribution at different times and experiments:

% Experiments to plot:
Exps     = [4 5 9];
Nexp_aux = numel(Exps);

% Times to plot (times of intersection between experiments):
ind45i = find(CFUST_ave_data(:,5) > CFUST_ave_data(:,4), 1);
ind45f = find(CFUST_ave_data(:,5) > CFUST_ave_data(:,4), 1, 'last');
ind49i = find(CFUST_ave_data(:,9) > CFUST_ave_data(:,4), 1);
ind49f = find(CFUST_ave_data(:,9) > CFUST_ave_data(:,4), 1, 'last');
ind59i = find(CFUST_ave_data(:,9) > CFUST_ave_data(:,5), 1);
ind59f = find(CFUST_ave_data(:,9) > CFUST_ave_data(:,5), 1, 'last');

tind_pdfAMR = sort([ind49i ind49f ind45i ind45f ind59i ind59f]);

nt_pdfAMR   = numel(tind_pdfAMR);


% ----------------------------------------------------------------------- %
% Options for area plot:

basestep = 0.15;                                                           % Set step between baselines;
baseline = zeros(nr, 1) + basestep*(nt_pdfAMR - 1);                        % Set base lines for area plots;
rtext    = rr([6 8 10 12 14 16]);                                          % x position of text in figure;
transp   = [1.0 0.6 0.6];                                                  % Transparency for the area plots;

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
    
            coord_up      = [rr, uplim];
            coord_low     = [rr, lowlim];
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Plot of deterministic VS stochastic heteroresistance distribution:

% ----------------------------------------------------------------------- %
% Calculate deterministic average:

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
rtext{1}    = rr([6 8 10 12 14 16]);
rtext{2}    = rr([6 8 10 12 14 16]);

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
    
            coord_up      = [rr, uplim];
            coord_low     = [rr, lowlim];
    
            coord_combine = [coord_up;flipud(coord_low)];
            
            subplot(1, 3, plot_count + 1)
            
            hold on
            
            fill(coord_combine(:,1),coord_combine(:,2), cc(iexp, :), 'EdgeColor', cccaux(plot_count,:), 'FaceAlpha', 1.0)
            
            % ----------------------------------------------- %
            % Plot model distribution:      
            uplim         = CFUS_mod(indt, 1:nr, iexp).'/CFUST_mod(indt, iexp) + lowlim;
            
            coord_up      = [rr, uplim];
    
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


print('-r720', 'PanelAMRlevelHD', '-dpng')





