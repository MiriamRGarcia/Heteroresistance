%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate confidence intervals using Fisher Information Matrix 
% for iid normal heterocedastic measurement noise using CFUS/mL as 
% the measured variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
clc

% ----------------------------------------------------------------------- %
% (1) Experimental setting:

% Confidence level:
confLevel  = 0.95;         

% Exact parameter values:
mugS    = 0.63;
mugR    = 0.36;
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

% Additional parameters for heterocedastic variance:
avar    = 1;
bvar    = 2;

% Array with parameter values:
par     = [mugS;mugR;alph_g;mukmaxS;bet;alphk;EC50k;Hk;xiSR;kxi;X0;lamb_IC;bvar];
nom_par = par;

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
rr  = linspace(ra, rb, nr)';
RR  = repmat(rr, 1, nr) - repmat(rr.', nr, 1);                              % Auxiliary matrix to solve ODEs,
RR  = RR - triu(RR) + tril(RR).';

% Parameters for input concentration:
MIC_S = EC50k*(mugS/(mukmaxS - mugS))^(1/Hk);
Cexp  = MIC_S*[0 1 2 4 8];
Nexp  = numel(Cexp);

% ----------------------------------------------------------------------- %
% Calculate total average counts xT and sensitivities:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

% Load optimal parameters:
Ntraj       = 3;
resultsname = sprintf('resFit_TD_Hetero_%dtraj_b2.mat', Ntraj);
load(resultsname, 'optSol', 'nom_par')

optPar      = optSol.*nom_par;
FIM_par     = optPar.';

% We scale with respect the true parameter values:
nom_par     = FIM_par; 

% Initial condition:
X0     = FIM_par(11);
lambIC = FIM_par(12);

f0 = exp(-lamb_IC*rr);
f0 = f0/sum(f0);
x0 = X0*f0;

% Obtain average and sensitivities:
[xT, sens_xT] = Sens_MultiExp(tmod, rr, Cexp, FIM_par, x0, ODEoptions);

% ----------------------------------------------------------------------- %
% (Optional) Plot sensitivities and normalised sensitivities if desired:
plot_sens = 0;

if plot_sens
    Cplot  = MIC_S*linspace(0, 8, 100).';
    Plot_Sens(tmod, Cplot, FIM_par, xT, sens_xT)
end

% ----------------------------------------------------------------------- %
% Calculate sensitivity matrix of model parameters:
SensMatrix = [];

np = numel(FIM_par) - 1;

for ip = 1:np
    
    % Sensitivity to parameter ip:
    aux = reshape(sens_xT(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    
    % Scaling by the nominal parameter (transformation par* = par/nom_par):
    aux = nom_par(ip)*aux;
    
    SensMatrix = [SensMatrix reshape(aux.', [], 1)];
end

% Rank of sensitivity matrix:
[~, SVS, ~]  = svd(SensMatrix);
rgSensMatrix = rank(SVS);

% Multi-exp normalised sensitivity matrix using the sampling times:      
% Important! this is the normalised sensitivity with respect parameters
% before reescaling. 
normSensMatrix = [];

for ip = 1:np
    aux = reshape(sens_xT(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    aux = FIM_par(ip)*aux./xT(texp_ind, 1:Nexp);
    normSensMatrix = [normSensMatrix reshape(aux.',[],1)];
end


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

% Calculate covariance matrix:
load(resultsname, 'CFUS_data')

CFUS_data = reshape(CFUS_data.', [], 1); 
xT        = reshape(xT(texp_ind, 1:Nexp).', [], 1); 

bvar = FIM_par(13);
avar = sum((CFUS_data - xT).^2./(xT.^bvar))/(ntexp*Nexp);

CovMatrix   = diag(avar*xT.^bvar);

[FIM, FConfInt] = Fisher_CI(FIM_par, avar, nom_par, SensMatrix, normSensMatrix, CovMatrix, confLevel);