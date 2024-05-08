function Plot_Sens(tmod, Cexp, par, xT, sens_xT)

% Problem sizes:
np   = size(sens_xT, 2);
nt   = numel(tmod);
Nexp = numel(Cexp);

% Parameter names:
name_pars = {'$\mu_{gS}$','$\mu_{gR}$','$\alpha_g$','$\mu_{k,max}^S$','$\beta$',...
             '$\alpha_k$','$EC_{50k}$','$H_k$','$\xi_{SR}$','$k_{\xi}$','$X_0$','$\lambda_{0}$'};

% Plot sensitivities:
figure(1)

for ip = 1:np
    
    subplot(4, 3, ip)
    
    xlabel('$t (h)$','Interpreter','Latex')
    ylabel('$C (mg/L)$','Interpreter','Latex')
    
    sens_aux = reshape(sens_xT(1:nt, ip, 1:Nexp), nt, Nexp);
    
    contourf(tmod, Cexp, sens_aux.', 10, 'edgecolor', 'none')

    CB   = colorbar;
    lCB  = get(CB,'Limits');
    tCB  = linspace(lCB(1),lCB(2),5);
    
    set(CB,'Ticks',tCB)
%     
%     TLCB = arrayfun(@(x) sprintf('%.2f',x), tCB, 'un', 0);
%     
%     set(CB,'TickLabels',TLCB)

    title(name_pars{ip}, 'Interpreter', 'Latex', 'FontSize', 11)
    
    xlabel('$t (h)$','Interpreter','Latex')
    ylabel('$C (mg/L)$','Interpreter','Latex')
    
    box off
    
end   

figure(2)

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



end