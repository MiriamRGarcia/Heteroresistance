%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate confidence intervals using Fisher Information Matrix for iid 
% normal measurement noise and known variance + deterministic data
% using CFUS/mL and the variances as measured variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

format long

% ----------------------------------------------------------------------- %
% Experimental setting:

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

% Array with parameter values:
par     = [mugS;mugR;alph_g;mukmaxS;bet;alphk;EC50k;Hk;xiSR;kxi;X0;lamb_IC];
nom_par = par;
np      = numel(par);

% Time discretisation:
t0       = 0; 
tf       = 48;       
ht       = 1e-2;             
tmod     = t0:ht:tf;        
texp     = [2 4 6 8 12 16 20 24 36 48];
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
MIC_S = 0.2;
Cexp  = MIC_S*[0 1 2 4 8];
Nexp  = numel(Cexp);

% Initial condition:
f0 = exp(-lamb_IC*rr);
f0 = f0/sum(f0);
x0 = X0*f0;

% ----------------------------------------------------------------------- %
% Calculate average and sensitivities:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

[xT, varxT, sens_xT, sens_varxT] = Sens_MultiExp(tmod, rr, Cexp, par, x0, ODEoptions);

% ----------------------------------------------------------------------- %
% (Optional) Plot sensitivities and normalised sensitivities if desired:
figs = 1;

if figs
    Plot_Sens(tmod, Cexp, par, xT, sens_xT, varxT, sens_varxT)
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

% Multi-exp normalised sensitivity matrix using the sampling times:       
scSensMatrix = [];
for ip = 1:np
    aux = reshape(sens_xT(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    aux = par(ip)*aux./xT(texp_ind, 1:Nexp);
    scSensMatrix = [scSensMatrix reshape(aux.',[],1)];
end
eig_Sens   = eig(SensMatrix.'*SensMatrix);
rSens      = rank(SensMatrix.'*SensMatrix);
eig_scSens = eig(scSensMatrix.'*scSensMatrix);
rscSens    = rank(scSensMatrix.'*scSensMatrix);

% Matrix with scaled sensitivities of variance:
SensMatrixVar   = [];
scSensMatrixVar = [];
for ip = 1:np
    aux = reshape(sens_varxT(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    % With scaling:
    SensMatrixVar   = [SensMatrixVar reshape(nom_par(ip)*aux.',[],1)];
    
    % Normalised sensitivity:
    aux = par(ip)*aux./varxT(texp_ind, 1:Nexp);
    scSensMatrixVar = [scSensMatrixVar reshape(aux.',[],1)];
end


% ----------------------------------------------------------------------- %
% Calculate correlation matrix:
cal_corr = 1;

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

% Calculate covariance matrix:
dgCovMatrix = reshape(varxT(texp_ind, 1:Nexp).', [], 1); 
CovMatrix   = diag(dgCovMatrix);

% Half-lenght of the confidence interval for each parameter:
FCI = Fisher_CI(Nexp, ntexp, par, nom_par, SensMatrix, SensMatrixVar, scSensMatrixVar, CovMatrix, cl);                    