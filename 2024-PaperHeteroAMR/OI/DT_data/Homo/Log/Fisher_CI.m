%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate confidence intervals using Fisher Information Matrix
% for normal iid error with homoscedastic measurement noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FIM, FConfInt] = Fisher_CI(Nexp, ntexp, par, SensMatrix, stDev, confLevel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np = numel(par);

% ----------------------------------------------------------%
% Scale the sensitivity matrix by the covariance matrix:
inv_stDev  = diag(ones(ntexp*Nexp, 1)/stDev);
SensMatrix = inv_stDev*SensMatrix;

% ----------------------------------------------------------%
% Singular value decomposition of Sensitivity Matrix:
[~, Sens_SingVal, Sens_VV] = svd(SensMatrix);

dSingVal      = diag(Sens_SingVal);
dSingVal_2    = dSingVal.^2;
SingVal_2     = diag(dSingVal_2);
InvdSingVal_2 = 1./dSingVal_2;
InvSingVal_2  = diag(InvdSingVal_2);

% ----------------------------------------------------------%
% Fisher Information and inverse:
FIM  = Sens_VV*SingVal_2*Sens_VV.';
IFIM = Sens_VV*InvSingVal_2*Sens_VV.';

% ----------------------------------------------------------%
% Calculate the tStudent value at the required confidence level:
alph  = 1 - confLevel;
pLo   = alph/2;
pUp   = 1 - alph/2;
ndf   = ntexp*Nexp - np;
tStud = tinv([pLo pUp], ndf);

% ----------------------------------------------------------%
% Calculate confidence interval:
pvar = diag(IFIM);

% With scaling:
%RFConfInt   = tStud(2)*sqrt(pvar)
%FConfInt    = nom_par.*RFConfInt

% Without scaling:
FConfInt   = tStud(2)*sqrt(pvar);


end