%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Fisher Information Matrix and confidence intervals using
% the Linear Noise Approximation for the approximate process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

%-------------------------------------------------------------------------%
% Experimental settings
%-------------------------------------------------------------------------%

% Parameter values:
mugS    = 0.63;
mugR    = 0.36;

alphg   = 2;

mukmaxS = 3.78;
bet     = 0.4;
alphk   = 3;

EC50k   = 1;
Hk      = 1;

xiSR    = 1e-6;
kxi     = log(1e2);

% Array of parameter values:
par = [mugS;mugR;alphg;mukmaxS;bet;alphk;EC50k;Hk;xiSR;kxi];

% Time discretisation:
t0       = 0; 
tf       = 48;       
ht       = 1e-3;             
tmod     = t0:ht:tf;        
%texp     = [2 4 6 8 12 16 20 24 36 48];
texp     = {[2 4 6 8 12 16 20 24 36 48],[2 4 6 8 12 16 20 24 36 48],[2 4 6 8 12 16 20 24 36 48],[2 4 6 8 12 16 20 24 36 48]};
tmod     = sort(unique([tmod texp{end}])).'; 
texp_ind = {find(ismember(tmod, texp{1})),find(ismember(tmod, texp{2})),find(ismember(tmod, texp{3})),find(ismember(tmod, texp{4}))} ;  
nt       = numel(tmod);    
%ntexp    = numel(texp);
ntexp   = [numel(texp{1}); numel(texp{2}); numel(texp{3}); numel(texp{4})];

% Discretisation of AMR level:
ra  = 0;
rb  = 1;
nr  = 50;
rr  = linspace(ra, rb, nr).';
RR  = repmat(rr, 1, nr) - repmat(rr.', nr, 1);                             % Auxiliary matrix to solve ODEs,
RR  = RR - triu(RR) + tril(RR).';

% Drug concentration in each experiment:
MIC_S = EC50k*(mugS/(mukmaxS - mugS))^(1/Hk);
Cexp = MIC_S*[0 0.1 0.2 0.3];
Nexp = numel(Cexp);

% Initial condition:
X0 = 1;
IC = 'Exp'; % 'S';

if strcmp(IC, 'S')
    % (Only S cells):
    f0  = eye(nr, 1);
    par = [par;X0];
else
    % Frecuency decaying exponentially:
    lamb_IC = 50;

    par = [par;X0;lamb_IC];

    f0  = exp(-lamb_IC*rr);
    f0  = f0/sum(f0);
end

x0 = X0*f0;

np = numel(par);

%-------------------------------------------------------------------------%
% Calculate average, variances and their sensitivities
%-------------------------------------------------------------------------%
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

[xT, varxT, sens_xT, sens_varxT, munT, sens_munT] = Sens_MV_MultiExp(tmod, rr, Cexp, par, x0, ODEoptions);

%-------------------------------------------------------------------------%
% Graphics (optional)
%-------------------------------------------------------------------------%
figs = 1;

if figs
    figure(1)
    hold on
    for iexp = 1:Nexp
        plot(tmod, log10(xT(1:nt, iexp)), 'DisplayName', strcat('C=',num2str(Cexp(iexp))))
    end
    legend show
    hold off
    xlabel('t')
    ylabel('CFUS/mL')
    
    figure(2)
    hold on
    for iexp = 1:Nexp
        plot(tmod, log10(varxT(1:nt, iexp)), 'DisplayName', strcat('C=',num2str(Cexp(iexp))))
    end
    legend show
    hold off
    xlabel('t')
    ylabel('\sigma^2')
    
    figure(3)
    
    sc_sens_xT = zeros(nt, np, Nexp);
    
    for ip = 1:np
        sens_aux = reshape(sens_xT(1:nt, ip, 1:Nexp), nt, Nexp);
        sc_sens_xT(1:nt, ip, 1:Nexp) = par(ip)*sens_aux./xT(1:nt, 1:Nexp);
        
        subplot(4, 3, ip)
        
        contourf(tmod, Cexp, reshape(sc_sens_xT(1:nt, ip, 1:Nexp), nt, Nexp).', 10, 'edgecolor', 'none')
        xlabel('$t$ (h)', 'Interpreter','Latex','FontSize', 11)
        ylabel('$C$ (mg/L)', 'Interpreter','Latex','FontSize', 11)

        CB   = colorbar;
        lCB  = get(CB,'Limits');
        tCB  = linspace(lCB(1),lCB(2),5);
        set(CB,'Ticks',tCB)
        TLCB = arrayfun(@(x) sprintf('%.2f',x), tCB, 'un', 0);
        set(CB,'TickLabels',TLCB)

        title(strcat('$S\theta$',num2str(ip)), 'Interpreter', 'Latex', 'FontSize', 11)
        xlim([0 48])
        box off
    end   
    
    figure(4)
    
    sc_sens_varxT = zeros(nt, np, Nexp);
    
    for ip = 1:np
        sens_aux = reshape(sens_varxT(1:nt, ip, 1:Nexp), nt, Nexp);
        sc_sens_varxT(1:nt, ip, 1:Nexp) = par(ip)*sens_aux./varxT(1:nt, 1:Nexp);
        
        subplot(4, 3, ip)
        
        contourf(tmod, Cexp, reshape(sc_sens_xT(1:nt, ip, 1:Nexp), nt, Nexp).', 10, 'edgecolor', 'none')
        xlabel('$t$ (h)', 'Interpreter','Latex','FontSize', 11)
        ylabel('$C$ (mg/L)', 'Interpreter','Latex','FontSize', 11)

        CB   = colorbar;
        lCB  = get(CB,'Limits');
        tCB  = linspace(lCB(1),lCB(2),5);
        set(CB,'Ticks',tCB)
        TLCB = arrayfun(@(x) sprintf('%.2f',x), tCB, 'un', 0);
        set(CB,'TickLabels',TLCB)

        title(strcat('$S\sigma^2\theta$',num2str(ip)), 'Interpreter', 'Latex', 'FontSize', 11)
        xlim([0 48])
        box off
    end         

end


%-------------------------------------------------------------------------%
% Calculate covariance matrix and sensitivities of covariances
%-------------------------------------------------------------------------%
Cov             = [];
inv_Cov         = [];
FI_xT           = [];

sensCov_mugS    = [];
sensCov_mugR    = [];
sensCov_alphg   = [];
sensCov_mukmaxS = [];
sensCov_bet     = [];
sensCov_alphk   = [];
sensCov_EC50k   = [];
sensCov_Hk      = [];
sensCov_xiSR    = [];
sensCov_kxi     = [];
sensCov_X0      = [];
sensCov_lamb_IC = [];

for iexp = 1:Nexp
    
    texp_ind_curr       = texp_ind{iexp};
    ntexp_curr          = ntexp(iexp);
    
    FI_xT               = [FI_xT;
                           xT(texp_ind_curr, iexp)];
                       
    aux_Cov             = diag(varxT(texp_ind_curr, iexp));               % Initialice with variances in the diagonal,
 
    aux_sensCov_mugS    = diag(sens_varxT(texp_ind_curr, 1, iexp));
    aux_sensCov_mugR    = diag(sens_varxT(texp_ind_curr, 2, iexp));
    aux_sensCov_alphg   = diag(sens_varxT(texp_ind_curr, 3, iexp));
    aux_sensCov_mukmaxS = diag(sens_varxT(texp_ind_curr, 4, iexp));
    aux_sensCov_bet     = diag(sens_varxT(texp_ind_curr, 5, iexp));
    aux_sensCov_alphk   = diag(sens_varxT(texp_ind_curr, 6, iexp));
    aux_sensCov_EC50k   = diag(sens_varxT(texp_ind_curr, 7, iexp));
    aux_sensCov_Hk      = diag(sens_varxT(texp_ind_curr, 8, iexp));
    aux_sensCov_xiSR    = diag(sens_varxT(texp_ind_curr, 9, iexp));
    aux_sensCov_kxi     = diag(sens_varxT(texp_ind_curr, 10, iexp));
    aux_sensCov_X0      = diag(sens_varxT(texp_ind_curr, 11, iexp));
    aux_sensCov_lamb_IC = diag(sens_varxT(texp_ind_curr, 12, iexp));
    
    for it = 1:(ntexp_curr - 1)
        z0         = [1;
                      zeros(np, 1)];
        taux       = tmod(texp_ind_curr(it):nt);
        texp_aux   = find(ismember(taux, texp{iexp}));
        
        [~, CVout] = ode15s(@(t, z) Odes_SensCV(t, z, tmod, munT(1:nt, iexp), sens_munT(1:nt, 1:np, iexp)), taux, z0, ODEoptions);
        
        cv      = CVout(texp_aux(2:end), 1);
        sens_cv = CVout(texp_aux(2:end), 2:end);
        
        aux_Cov(it, (it + 1):ntexp_curr) = varxT(texp_ind_curr(it))*cv.';
        
        aux_sensCov_mugS(it, (it + 1):ntexp_curr)    = sens_varxT(texp_ind_curr(it), 1, iexp)*cv + varxT(texp_ind_curr(it))*sens_cv(:,1);
        aux_sensCov_mugR(it, (it + 1):ntexp_curr)    = sens_varxT(texp_ind_curr(it), 2, iexp)*cv + varxT(texp_ind_curr(it))*sens_cv(:,2);
        aux_sensCov_alphg(it, (it + 1):ntexp_curr)   = sens_varxT(texp_ind_curr(it), 3, iexp)*cv + varxT(texp_ind_curr(it))*sens_cv(:,3);
        aux_sensCov_mukmaxS(it, (it + 1):ntexp_curr) = sens_varxT(texp_ind_curr(it), 4, iexp)*cv + varxT(texp_ind_curr(it))*sens_cv(:,4);
        aux_sensCov_bet(it, (it + 1):ntexp_curr)     = sens_varxT(texp_ind_curr(it), 5, iexp)*cv + varxT(texp_ind_curr(it))*sens_cv(:,5);
        aux_sensCov_alphk(it, (it + 1):ntexp_curr)   = sens_varxT(texp_ind_curr(it), 6, iexp)*cv + varxT(texp_ind_curr(it))*sens_cv(:,6);
        aux_sensCov_EC50k(it, (it + 1):ntexp_curr)   = sens_varxT(texp_ind_curr(it), 7, iexp)*cv + varxT(texp_ind_curr(it))*sens_cv(:,7);
        aux_sensCov_Hk(it, (it + 1):ntexp_curr)      = sens_varxT(texp_ind_curr(it), 8, iexp)*cv + varxT(texp_ind_curr(it))*sens_cv(:,8);
        aux_sensCov_xiSR(it, (it + 1):ntexp_curr)    = sens_varxT(texp_ind_curr(it), 9, iexp)*cv + varxT(texp_ind_curr(it))*sens_cv(:,9);
        aux_sensCov_kxi(it, (it + 1):ntexp_curr)     = sens_varxT(texp_ind_curr(it), 10, iexp)*cv + varxT(texp_ind_curr(it))*sens_cv(:,10);
        aux_sensCov_X0(it, (it + 1):ntexp_curr)      = sens_varxT(texp_ind_curr(it), 11, iexp)*cv + varxT(texp_ind_curr(it))*sens_cv(:,11);
        aux_sensCov_lamb_IC(it, (it + 1):ntexp_curr) = sens_varxT(texp_ind_curr(it), 12, iexp)*cv + varxT(texp_ind_curr(it))*sens_cv(:,12);
        
    end
    
    aux_Cov     = aux_Cov + triu(aux_Cov, 1).';
    aux_inv_Cov = inv(aux_Cov);
    
    aux_sensCov_mugS    = aux_sensCov_mugS + triu(aux_sensCov_mugS, 1).';
    aux_sensCov_mugR    = aux_sensCov_mugR + triu(aux_sensCov_mugR, 1).';
    aux_sensCov_alphg   = aux_sensCov_alphg + triu(aux_sensCov_alphg, 1).';
    aux_sensCov_mukmaxS = aux_sensCov_mukmaxS + triu(aux_sensCov_mukmaxS, 1).';
    aux_sensCov_bet     = aux_sensCov_bet + triu(aux_sensCov_bet, 1).';
    aux_sensCov_alphk   = aux_sensCov_alphk + triu(aux_sensCov_alphk, 1).';
    aux_sensCov_EC50k   = aux_sensCov_EC50k + triu(aux_sensCov_EC50k, 1).';
    aux_sensCov_Hk      = aux_sensCov_Hk + triu(aux_sensCov_Hk, 1).';
    aux_sensCov_xiSR    = aux_sensCov_xiSR + triu(aux_sensCov_xiSR, 1).';
    aux_sensCov_kxi     = aux_sensCov_kxi + triu(aux_sensCov_kxi, 1).';
    aux_sensCov_X0      = aux_sensCov_X0 + triu(aux_sensCov_X0, 1).';
    aux_sensCov_lamb_IC = aux_sensCov_lamb_IC + triu(aux_sensCov_lamb_IC, 1).';
    
    Cov     = blkdiag(Cov, aux_Cov);
    inv_Cov = blkdiag(inv_Cov, aux_inv_Cov);
    
    sensCov_mugS    = blkdiag(sensCov_mugS, aux_sensCov_mugS);
    sensCov_mugR    = blkdiag(sensCov_mugR, aux_sensCov_mugR);
    sensCov_alphg   = blkdiag(sensCov_alphg, aux_sensCov_alphg);
    sensCov_mukmaxS = blkdiag(sensCov_mukmaxS, aux_sensCov_mukmaxS);
    sensCov_bet     = blkdiag(sensCov_bet, aux_sensCov_bet);
    sensCov_alphk   = blkdiag(sensCov_alphk, aux_sensCov_alphk);
    sensCov_EC50k   = blkdiag(sensCov_EC50k, aux_sensCov_EC50k);
    sensCov_Hk      = blkdiag(sensCov_Hk, aux_sensCov_Hk);
    sensCov_xiSR    = blkdiag(sensCov_xiSR, aux_sensCov_xiSR);
    sensCov_kxi     = blkdiag(sensCov_kxi, aux_sensCov_kxi);
    sensCov_X0      = blkdiag(sensCov_X0, aux_sensCov_X0);
    sensCov_lamb_IC = blkdiag(sensCov_lamb_IC, aux_sensCov_lamb_IC);
end

sensCov(:, :, 1) = sensCov_mugS;
sensCov(:, :, 2) = sensCov_mugR;
sensCov(:, :, 3) = sensCov_alphg;
sensCov(:, :, 4) = sensCov_mukmaxS;
sensCov(:, :, 5) = sensCov_bet;
sensCov(:, :, 6) = sensCov_alphk;
sensCov(:, :, 7) = sensCov_EC50k;
sensCov(:, :, 8) = sensCov_Hk;
sensCov(:, :, 9) = sensCov_xiSR;
sensCov(:, :, 10) = sensCov_kxi;
sensCov(:, :, 11) = sensCov_X0;
sensCov(:, :, 12) = sensCov_lamb_IC;

%-------------------------------------------------------------------------%
% Calculate confidence intervals using Fisher Information Matrix:
cl  = 0.9;                                                                % Confidence level,

%sens_xT = sens_xT(texp_ind, 1:np, 1:Nexp);
sens_xT_FI = [];

for iexp = 1:Nexp
    sens_aux   = reshape(sens_xT(texp_ind{iexp}, 1:np, iexp), ntexp(iexp), np);
    sens_xT_FI = [sens_xT_FI;
                  sens_aux];
end

FCI = Fisher_CI(sens_xT_FI, Cov, inv_Cov, sensCov, cl, Nexp);    

