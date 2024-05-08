function Plot_Sens(tmod, Cexp, par, xT, sens_xT, varxT, sens_varxT)

% Problem sizes:
np   = size(sens_xT, 2);
nt   = numel(tmod);
Nexp = numel(Cexp);

% Parameter names:
name_pars = {'$\mu_{gS}$','$\mu_{gR}$','$\alpha_g$','$\mu_{k,max}^S$','$\beta$',...
             '$\alpha_k$','$EC_{50k}$','$H_k$','$\xi_{SR}$','$k_{\xi}$','$X_0$','$\lambda_{0}$'};

% ----------------------------------------------------------------------- %
% Plot sensitivities of the average:
ifig = 1;
% figure(ifig)
% 
% for ip = 1:np
%     
%     subplot(4, 3, ip)
%     
%     xlabel('$t (h)$','Interpreter','Latex')
%     ylabel('$C (mg/L)$','Interpreter','Latex')
%     
%     sens_aux = reshape(sens_xT(1:nt, ip, 1:Nexp), nt, Nexp);
%     
%     contourf(tmod, Cexp, sens_aux.', 10, 'edgecolor', 'none')
% 
%     CB   = colorbar;
%     lCB  = get(CB,'Limits');
%     tCB  = linspace(lCB(1),lCB(2),5);
%     
%     set(CB,'Ticks',tCB)
% %     
% %     TLCB = arrayfun(@(x) sprintf('%.2f',x), tCB, 'un', 0);
% %     
% %     set(CB,'TickLabels',TLCB)
% 
%     title(name_pars{ip}, 'Interpreter', 'Latex', 'FontSize', 11)
%     
%     xlabel('$t (h)$','Interpreter','Latex')
%     ylabel('$C (mg/L)$','Interpreter','Latex')
%     
%     box off
%     
% end  
% ifig = ifig + 1;

% ----------------------------------------------------------------------- %
% Plot scaled sensitivities of the average:
figure(ifig)

for ip = 1:np
    
    subplot(4, 3, ip)
    
    sens_aux    = reshape(sens_xT(1:nt, ip, 1:Nexp), nt, Nexp);
    sc_sens_aux = par(ip)*sens_aux./xT;
    
    contourf(tmod, Cexp, sc_sens_aux.', 10, 'edgecolor', 'none')

    CB   = colorbar;
    lCB  = get(CB,'Limits');
    tCB  = linspace(lCB(1),lCB(2),5);
    
    set(CB,'Ticks',tCB)
    
%     TLCB = arrayfun(@(x) sprintf('%.2f',x), tCB, 'un', 0);
%     
%     set(CB,'TickLabels',TLCB)

    title(name_pars{ip}, 'Interpreter', 'Latex', 'FontSize', 11)
    
    xlabel('$t (h)$','Interpreter','Latex')
    ylabel('$C (mg/L)$','Interpreter','Latex')
    
    box off
end

ifig = ifig + 1;

% ----------------------------------------------------------------------- %
% Plot scaled sensitivities of the standard deviation:
figure(ifig)
    
sc_sens_varxT = zeros(nt, np, Nexp);
    
for ip = 1:np
    sens_aux = reshape(sens_varxT(1:nt, ip, 1:Nexp), nt, Nexp);
    sc_sens_varxT(1:nt, ip, 1:Nexp) = par(ip)*sens_aux./varxT(1:nt, 1:Nexp);
        
    subplot(4, 3, ip)
        
    contourf(tmod, Cexp, reshape(sc_sens_varxT(1:nt, ip, 1:Nexp), nt, Nexp).', 10, 'edgecolor', 'none')
    
    xlabel('$t$ (h)', 'Interpreter','Latex','FontSize', 11)
    ylabel('$C$ (mg/L)', 'Interpreter','Latex','FontSize', 11)

    CB   = colorbar;
    lCB  = get(CB,'Limits');
    tCB  = linspace(lCB(1),lCB(2),5);
    set(CB,'Ticks',tCB)
%     TLCB = arrayfun(@(x) sprintf('%.2f',x), tCB, 'un', 0);
%     set(CB,'TickLabels',TLCB)

    title(name_pars{ip}, 'Interpreter', 'Latex', 'FontSize', 11)

    box off       
end



end