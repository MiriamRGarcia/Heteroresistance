%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main script to calculate Fisher Confidence Intervals  
% for iid normal measurement homoscedastic noise with log10(CFUS/mL) data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
clc

% ----------------------------------------------------------------------- %
% (1) Experimental setup:

% Confidence level:
confLevel  = 0.95;         

% Load optimal fit parameters:
Ntraj       = 3;

resultsname = sprintf('resFit_TD_%dtraj_log.mat', Ntraj);
load(resultsname, 'optPar')

optPar      = optPar(end, :).';
FIM_par     = optPar;                                                      

% Exact parameter values:
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

X0      = 1e6;
lambIC  = 50;

% Standard deviation of measurement noise:
stDev_exact = 0.5;

% Exact parameter array:
par     = [mugS;mugR;alphg;mukmaxS;bet;alphk;EC50k;Hk;xiSR;kxi;X0;lambIC];
nom_par = [0.1;0.1;1;1;0.1;1;1;1;1e-6;1;1e6;10];
np      = numel(par);

% Time discretisation:
t0        = 0; 
tf        = 48;       
ht        = 1e-3;             
tmod      = t0:ht:tf;        
texp      = [2 4 6 8 10 12 16 20 24 36 48];
tmod      = sort(unique([tmod texp])).'; 
texp_ind  = find(ismember(tmod, texp));  
nt        = numel(tmod);    
ntexp     = numel(texp);

% Discretisation of AMR level:
ra  = 0;
rb  = 1;
nr  = 50;
rr  = linspace(ra, rb, nr)';
RR  = repmat(rr, 1, nr) - repmat(rr.', nr, 1);                             % Auxiliary matrix to solve ODEs,
RR  = RR - triu(RR) + tril(RR).';

% Drug concentration in each experiment:
MIC_S = EC50k*(mugS/(mukmaxS - mugS))^(1/Hk);
Cexp  = MIC_S*[0.0 1.0 2.0 4.0 8.0].';
Nexp  = numel(Cexp);

FIM_par = par;
% ----------------------------------------------------------------------- %
% (2) Calculate total average counts xT and sensitivities:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

X0      = FIM_par(11);
lambIC  = FIM_par(12);

% Initial condition:
f0  = exp(-lambIC*rr);
f0  = f0/sum(f0);
x0  = X0*f0;

% Function to calculate sensitivities of xT:
[xT, sens_xT] = Sens_MultiExp(tmod, rr, Cexp, FIM_par, x0, ODEoptions);

% ----------------------------------------------------------------------- %
% (3) Calculat Log-10 measures of the average xT and sensitivities:

% Log10 of total average count:
yT = log10(xT);

% Calculate sensitivities of yT to parameters:
sens_yT = zeros(nt, np, Nexp);

for ip = 1:np
    aux = reshape(sens_xT(1:nt, ip, 1:Nexp), nt, Nexp);
    aux = (aux./xT)/log(10);
    sens_yT(1:nt, ip, 1:Nexp) = aux;
end

% ----------------------------------------------------------------------- %
% (4) (Optional) Plot normalised sensitivities if desired:
plot_sens = 0;

if plot_sens
    Cplot  = MIC_S*linspace(0, 8, 50).';
    Plot_Sens(tmod, Cexp, FIM_par, yT, sens_yT)
    return
end

% ----------------------------------------------------------------------- %
% (5) Calculate sensitivity matrix:
SensMatrix = [];
for ip = 1:np
    
    % Obtain sensitivities at the sampling times:
    aux = reshape(sens_yT(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    
    % Scaling by the nominal parameter (transformation par* = par/nom_par):
    %aux = nom_par(ip)*aux;
    
    % Scaling by the log ofparameter (transformation par* = log10(par)):
    %aux = FIM_par(ip)*log(10)*aux;
    
    SensMatrix = [SensMatrix reshape(aux.', [], 1)];
end

% Rank of sensitivity matrix:
[~, SVS, VS] = svd(SensMatrix);
rgSensMatrix = rank(SVS);

% Multi-exp normalised sensitivity matrix using the sampling times:       
normSensMatrix = [];
for ip = 1:np
    aux = reshape(sens_yT(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    aux = FIM_par(ip)*aux./yT(texp_ind, 1:Nexp);
    normSensMatrix = [normSensMatrix reshape(aux.',[],1)];
end
return
% ----------------------------------------------------------------------- %
% (6) (Optional) Calculate correlation matrix:
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
% (7) Calculate confidence intervals from Fisher Information Matrix:

% Obtain estimate for the homoscedastic variance:
load(resultsname, 'logCFUS_data')

logCFUS_data = reshape(logCFUS_data, [], 1);
yT           = reshape(yT(texp_ind, 1:Nexp), [], 1);

optVar       = sum((logCFUS_data - yT).^2)/(ntexp*Nexp);
optDev       = sqrt(optVar);

[FIM, FConfInt] = Fisher_CI(Nexp, ntexp, FIM_par, SensMatrix, optDev, confLevel);  




