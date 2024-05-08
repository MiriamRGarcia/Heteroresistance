%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot_Sens(tmod, Cexp, par, xT, sens_xT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% General setup:

% Problem sizes:
np   = size(sens_xT, 2);
nt   = numel(tmod);
Nexp = numel(Cexp);

% Parameter names:
name_pars = {'$\mu_{bS}$','$\mu_{bR}$','$\alpha_b$','$\mu_{dmax,S}$','$\alpha_d$',...
             '$\beta_d$','$EC_{50d}$','$H_d$','$\xi_{SR}$','$k_{\xi}$','$N_{T0}$','$\lambda_{T0}$'};

% Changue order of alphk and betk:
paraux5   = par(5);

par(5)    = par(6);
par(6)    = paraux5;

aux_alphk = sens_xT(1:nt, 6, 1:Nexp);
aux_betk  = sens_xT(1:nt, 5, 1:Nexp);

sens_xT(1:nt, 5, 1:Nexp) = aux_alphk;
sens_xT(1:nt, 6, 1:Nexp) = aux_betk;


% ----------------------------------------------------------------------- %
% Plot normalised sensitivities:

fig = figure;

set(gcf,'color','w');

fig.Position = [554 152 1128 783];

Cticks = linspace(Cexp(1), Cexp(Nexp), 4);

for ip = 1:np
    
    subplot(4, 3, ip)
   
    % Calculate normalised sensitivities:
    sens_aux    = reshape(sens_xT(1:nt, ip, 1:Nexp), nt, Nexp);
    sc_sens_aux = par(ip)*sens_aux./xT;
    
    sc_sens_aux = abs(sc_sens_aux);
    
    % Intermediate step to create common colorbar:
    % smin = min(smin, min(min(sc_sens_aux)))
    % smax = max(smax, max(max(sc_sens_aux)))
    

    % Contour plot:
    contourf(tmod, Cexp, sc_sens_aux.', 100, 'edgecolor', 'none')
    
    
    % Set scale for parameter N0:
    if ip == 11
        caxis([0 max(max(sc_sens_aux))]);
    end
    
    % Style:
    
    %colormap(copper)
    %colormap(hot)
    %colormap(summer)
    
    % Colorbar of each plot:
    CB   = colorbar;
    lCB  = get(CB,'Limits');
    tCB  = linspace(lCB(1), lCB(2), 4);
    
    CB.FontSize = 11;
    
    set(CB,'Ticks',tCB,'TickLabelInterpreter', 'Latex')
    
    TLCB = arrayfun(@(x) sprintf('%.2f',x), tCB, 'un', 0);
    
    set(CB,'TickLabels',TLCB)
    
    
    % Axes limits:    
    xlim([tmod(1) tmod(end)])
    ylim([Cexp(1) Cexp(end)])
    
    % Axis properties:
    set(gca,'XTick',[0 10 20 30 40 48], 'YTick', Cticks, 'TickLabelInterpreter', 'Latex', 'FontSize', 12)
    
    % Title:
    title(name_pars{ip}, 'Interpreter', 'Latex', 'FontSize', 16)
    
    % Set ticks on y axis:
    ax = gca;
    ax.YAxis.TickLabelFormat = '%.2f';
    
    box off
end

% Only one xy axes for the figure:
han = axes(fig,'visible','off'); 

han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';


Xlim = xlim;           
Xlb  = mean(Xlim); 
xlabel(han,{'','Time (h)'}, 'Interpreter','Latex','FontSize', 16, 'Position',[Xlb-0.03 -0.07],'HorizontalAlignment','center')

ylabel(han,{'$C$ (mg/L)',''}, 'Interpreter','Latex','FontSize', 16);

% Export figure:
%export_fig('Panel_Sens_ave','-jpg','-transparent')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two colormaps: 

% col1 = cmocean('amp',5); % 6 shades of red
% col2 = parula(5);       % 95 shades of parula
% % Concatenate colormaps: 
% col = cat(1,col1,col2); 
% colormap(col)


% cmap       = colormap  %current map
% cmap(51,:) = [1 1 1];  %make first color white
% colormap(cmap);        %put it into effect

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common colorbar:
% CB   = colorbar;
% tCB  = linspace(smin, smax, 10);
%     
% set(CB,'Ticks',tCB)
%     
% TLCB = arrayfun(@(x) sprintf('%.2f',x), tCB, 'un', 0);
%     
% set(CB,'TickLabels',TLCB)
% 
% hp4 = get(subplot(4, 3, 12),'Position');
% colorbar('Position', [hp4(1) + hp4(3) + 0.01  hp4(2)  0.02  hp4(1) + 0.05])

end