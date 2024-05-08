%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate confidence intervals using Fisher Information Matrix
% for normal iid error with heterocedatic noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FIM, FConfInt] = Fisher_CI(par, avar, nom_par, SensMatrix, normSensMatrix, CovMatrix, confLevel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of parameters:
np   = numel(par);

% Obtain parameter of variance:
bvar = par(np);

ndata = size(CovMatrix, 1);

% ---------------------------------------------------------- %
% Calculate the term in the sensitivity matrix:

% Scale the sensitivity matrix by the measurement standard deviations:
dgCov  = diag(CovMatrix);
dgSD   = sqrt(dgCov);
inv_SD = diag(1./dgSD);
SS     = inv_SD*SensMatrix;


% Ampliate with derivatives of the average to bvar:
SS     = [SS zeros(ndata, 1)];

% Singular value decomposition of SS:
[~, SS_SV, SS_V] = svd(SS);

dSS_SV = diag(SS_SV);     % Array with singular values of SS,
SS_SV2 = diag(dSS_SV.^2); % Diagonal matrix with square of singular values,

% First term of Fisher Information Matrix:
FIM    = SS_V*SS_SV2*SS_V.';

% ----------------------------------------------------------%
% Calculate the term in covariances:
% Note: this is true taking the rescaled parameters p^* = p/ptrue

%scSS  = [bvar*normSensMatrix ones(ntexp*Nexp, 1) bvar*log(dgCov)];
aux   = (dgCov/avar).^(1/bvar);
scSS  = [bvar*normSensMatrix bvar*log(aux)];

% Singular value decomposition of the scaled sensitivity matrix:
[~, scSS_SV, scSS_V] = svd(scSS);

dscSS_SV = diag(scSS_SV);     % Array with singular values of scSS,
scSS_SV2 = diag(dscSS_SV.^2); % Diagonal matrix with square of singular values,

TrCov = scSS_V*scSS_SV2*scSS_V.';

% Term of the FIM in the traces of covariance matrix sensitivities:
FIM = FIM + 0.5*TrCov;

% ----------------------------------------------------------%
% Calculate the inverse of FIM:

% Cholesky factorisation of FIM:
Chol_FIM = chol(FIM);

% Inverse of the FIM from svd of RR_FIM:
[~, Chol_FIM_SV, Chol_FIM_V] = svd(Chol_FIM);

dChol_FIM_SV  = diag(Chol_FIM_SV);           % Array with singular values of RR_FIM,
Chol_FIM_ISV2 = diag(1./(dChol_FIM_SV.^2));  % Diagonal matrix with inverse of squares of singular values,

IFIM = Chol_FIM_V*Chol_FIM_ISV2*Chol_FIM_V.';

% ----------------------------------------------------------%
% Calculate the tStudent value at the required confidence level:
alph = 1 - confLevel;
pLo  = alph/2;
pUp  = 1 - alph/2;
ndf  = ndata - np;
tStud = tinv([pLo pUp], ndf);

% ----------------------------------------------------------%
% Calculate confidence interval:
pvar = diag(IFIM);

% With scaling:
RFConfInt   = tStud(2)*sqrt(pvar);
FConfInt    = nom_par.*RFConfInt;


end