%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate confidence intervals using Fisher Information Matrix  
% for process noise using (CFUS/mL) as the measured variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
clc

% ----------------------------------------------------------------------- %
% (1) Experimental setting:

% Confidence level:
confLevel  = 0.95;         

% Load optimal parameters:
Ntraj   = 3;
resultsname = sprintf('resFit_TP_%dtraj.mat', Ntraj);
load(resultsname, 'optSol','nom_par', 'CFUS_data')
optPar  = optSol.*nom_par;
FIM_par = optPar.';
nom_par = FIM_par;

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

% Exact parameter array:
par     = [mugS;mugR;alphg;mukmaxS;bet;alphk;EC50k;Hk;xiSR;kxi;X0;lambIC];
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
% (3) (Optional) Plot normalised sensitivities if desired:
plot_sens = 0;

if plot_sens
    Cplot  = MIC_S*linspace(0, 8, 100).';
    Plot_Sens(tmod, Cplot, FIM_par, xT, sens_xT)
    return
end

% ----------------------------------------------------------------------- %
% (4) Calculate sensitivity matrix:
SensMatrix = [];

ind0 = find(CFUS_data == 0); % Remove time points with no data,

for ip = 1:np
    
    % Obtain sensitivities at the sampling times:
    aux = reshape(sens_xT(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    
    aux(ind0) = [];

    % Scaling by the nominal parameter (transformation par* = par/nom_par):
    aux = nom_par(ip)*aux;
    
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
    aux = reshape(sens_xT(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    aux = FIM_par(ip)*aux./xT(texp_ind, 1:Nexp);
    normSensMatrix = [normSensMatrix reshape(aux.',[],1)];
end

% ----------------------------------------------------------------------- %
% (5) (Optional) Calculate correlation matrix:
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
% (6) Calculate covariance matrix:
load(resultsname, 'Var_data')
Var_data(Var_data == 1) = [];
Var_data = reshape(Var_data, [], 1);
CovMatrix  = diag(Var_data);


% ----------------------------------------------------------------------- %
% (7) Calculate confidence intervals from Fisher Information Matrix:
[FIM, FConfInt] = Fisher_CI(FIM_par, nom_par, SensMatrix, CovMatrix, confLevel);  




