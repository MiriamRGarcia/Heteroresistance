%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate confidence intervals using Fisher Information Matrix for iid 
% normal measurement noise and known variance + deterministic data
% using CFUS/mL as measured variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

format long

% ----------------------------------------------------------------------- %
% Experimental setting:

par_val = 'E'; %'E'                                                        % Use exact or estimate value to calculate FIM,

% Parameter values (E.Coli):
mugS    = 0.63;%log(2)*3; %[h-1]
mugR    = 0.36;%0.96*mugS;
alph_g  = 2;

mukmaxS = 3.78;
bet     = 0.4;
alphk   = 3;

EC50k   = 1;
Hk      = 1;

xiSR    = 1e-6;
xiM     = 1e-2;
kxi     = log(1e2);

X0      = 1e6;
lamb_IC = 50;

par     = [mugS;mugR;alph_g;mukmaxS;bet;alphk;EC50k;Hk;xiSR;kxi;X0;lamb_IC];
nom_par = [0.1;0.1;1;1;0.1;1;1;1;1e-6;1;1e6;10];
%nom_par = par;

if strcmp(par_val, 'L')   
    load('resFit_TD_3.mat','optSol', 'nom_par')
    
    par = (optSol.*nom_par).';
    
    mugS    = par(1);
    mugR    = par(2);
    alph_g  = par(3);

    mukmaxS = par(4);
    bet     = par(5);
    alphk   = par(6);

    EC50k   = par(7);
    Hk      = par(8);

    xiSR    = par(9);
    kxi     = par(10);

    X0      = par(11);
    lamb_IC = par(12);
    
    nom_par = nom_par.';
    %nom_par = par;
end
    
np = numel(par);

% Time discretisation:
t0       = 0; 
tf       = 48;       
ht       = 1e-3;             
tmod     = t0:ht:tf;        
texp     = [2 4 6 8 10 12 16 20 24 36 48];
tmod     = sort(unique([tmod texp])).'; 
texp_ind = find(ismember(tmod, texp));  
nt       = numel(tmod);    
ntexp    = numel(texp);

% Discretisation of the AMR level:
ra  = 0;
rb  = 1;
nr  = 50;
rr  = linspace(ra, rb, nr).';
RR  = repmat(rr, 1, nr) - repmat(rr.', nr, 1);                              % Auxiliary matrix to solve ODEs,
RR  = RR - triu(RR) + tril(RR).';

% Parameters for input concentration:
MIC_S = EC50k*(mugS/(mukmaxS - mugS))^(1/Hk);
Cexp  = MIC_S*[0 1 2 4 8];%linspace(0, 8, 100);
Nexp  = numel(Cexp);

% Initial condition:
f0 = exp(-lamb_IC*rr);
f0 = f0/sum(f0);
x0 = X0*f0;

% ----------------------------------------------------------------------- %
% Calculate average and sensitivities:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

[xT, sens_xT] = Sens_MultiExp(tmod, rr, Cexp, par, x0, ODEoptions);

% ----------------------------------------------------------------------- %
% (Optional) Plot sensitivities and normalised sensitivities if desired:
figs = 0;

if figs
    Plot_Sens(tmod, Cexp, par, xT, sens_xT)
    return
end

% ----------------------------------------------------------------------- %
% Calculate sensitivity matrix:
SensMatrix = [];

for ip = 1:np
    aux = reshape(sens_xT(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    
    % With scaling:
    aux = aux*nom_par(ip);
    SensMatrix = [SensMatrix reshape(aux.',[],1)];
end

[~, SVS, VS] = svd(SensMatrix);
rgS          = rank(SVS);



% Multi-exp normalised sensitivity matrix using the sampling times:       
scSensMatrix = [];

for ip = 1:np
    aux = reshape(sens_xT(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    aux = par(ip)*aux./xT(texp_ind, 1:Nexp);
    scSensMatrix = [scSensMatrix reshape(aux.',[],1)];
end

[~,SVscS, VscS] = svd(scSensMatrix);
rgscS           = rank(SVscS);

return
% ----------------------------------------------------------------------- %
% Calculate correlation matrix:
cal_corr = 0;

if cal_corr 
    par_corr = eye(np, np);
    for ip = 1:np
        for jp = 1:np
            zz1 = SensMatrix(:, ip);
            zz2 = SensMatrix(:, jp);
        
            sd1   = sqrt((zz1 - mean(zz1)).'*(zz1 - mean(zz1))/(ntexp*Nexp - 1));
            sd2   = sqrt((zz2 - mean(zz2)).'*(zz2 - mean(zz2))/(ntexp*Nexp - 1));
            cov12 = (zz1 - mean(zz1)).'*(zz2 - mean(zz2))/(ntexp*Nexp - 1);
        
            par_corr(ip, jp) = cov12/(sd1*sd2);
        end
    end
end

% ----------------------------------------------------------------------- %
% Calculate confidence intervals from Fisher Information Matrix:
cl  = 0.95;                                                                % Confidence level,
sd  = sqrt(1e4);                                                           % Constant standard deviation of measures,

FCI = Fisher_CI(Nexp, ntexp, par, nom_par, SensMatrix, sd, cl);            % Half-lenght of the confidence interval for each parameter,