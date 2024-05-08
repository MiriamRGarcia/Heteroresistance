%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate confidence intervals using Fisher Information Matrix
% for normal iid error with heterocedatic noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FCI = Fisher_CI(Nexp, ntexp, par, nom_par, SensMatrix, SensMatrixVar, scSensMatrixVar, CovMatrix, cl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of parameters:
np = numel(par);

% ---------------------------------------------------------- %
% Calculate the term in the sensitivity matrix:

% Scale the sensitivity matrix by the measurement standard deviations:
dgCov  = diag(CovMatrix);
dgSD   = sqrt(dgCov);
inv_SD = diag(1./dgSD);
SS     = inv_SD*SensMatrix;

% Singular value decomposition of SS:
[~, SS_SV, SS_V] = svd(SS);

dSS_SV = diag(SS_SV);                                                      % Array with singular values of SS,
SS_SV2 = diag(dSS_SV.^2);                                                  % Diagonal matrix with square of singular values,

% First term of Fisher Information Matrix:
FIM = SS_V*SS_SV2*SS_V.';

% ----------------------------------------------------------%
% Calculate the term in sensitivity of variances:

% [~, SVar_SV, SVar_V] = svd(SensMatrixVar);
% 
% dVar_SV  = diag(SVar_SV);                                                  % Array with singular values of scVar,
% SVar_SV2 = diag(dVar_SV.^2);                                               % Diagonal matrix with square of singular values,
% SensVar2 = SVar_V*SVar_SV2*SVar_V.';
% 
% FIM = FIM + SensVar2;
[~, SVar_SV, SVar_V] = svd(scSensMatrixVar);

dVar_SV  = diag(SVar_SV);                                                  % Array with singular values of scVar,
SVar_SV2 = diag(dVar_SV.^2);                                               % Diagonal matrix with square of singular values,
SensVar2 = SVar_V*SVar_SV2*SVar_V.';

FIM = FIM + SensVar2;

% ----------------------------------------------------------%
% Calculate the term in covariances:

% Singular value decomposition of the scaled sensitivity matrix for
% variances
% Note: this is true taking the rescaled parameters p^* = p/p_n
[~, scVar_SV, scVar_V] = svd(scSensMatrixVar);

dscVar_SV = diag(scVar_SV);                                                % Array with singular values of scVar,
scVar_SV2 = diag(dscVar_SV.^2);                                            % Diagonal matrix with square of singular values,

% Term in sensitivity matrix of variances as measured variable:
TrCov = scVar_V*scVar_SV2*scVar_V.'; 
FIM   = FIM + 0.5*TrCov;

% ----------------------------------------------------------%
% Calculate the inverse of FIM:

% Cholesky factorisation of FIM:
RR_FIM = chol(FIM);

% Inverse of the FIM from svd of RR_FIM:
[~, RR_FIM_SV, RR_FIM_V] = svd(RR_FIM);

dRR_FIM_SV  = diag(RR_FIM_SV);           % Array with singular values of RR_FIM,
RR_FIM_ISV2 = diag(1./(dRR_FIM_SV.^2));  % Diagonal matrix with inverse of squares of singular values,

IFIM = RR_FIM_V*RR_FIM_ISV2*RR_FIM_V.';

% Variances of parameters based on the inverse of FIM:
pvar = diag(IFIM);


% ----------------------------------------------------------%
% Calculate the tStudent value at the required confidence level:
alph = 1 - cl;
pLo  = alph/2;
pUp  = 1 - alph/2;
ndf  = ntexp*Nexp - np;
tStud = tinv([pLo pUp], ndf);

% ----------------------------------------------------------%
% Calculate confidence interval:

fprintf('\n >> El numero de condicion de la Fisher es:')
ncFIM  = cond(FIM)

fprintf('\n >> El rango de la Fisher es:')
rgFIM  = rank(RR_FIM_SV)

fprintf('\n >> El determinante de la Fisher es:')
detFIM = prod(dRR_FIM_SV)

fprintf('\n >> La mitad del ancho del IC (OJO, RELATIVO!!) es:')
FRCI   = tStud(2)*sqrt(pvar)

fprintf('\n >> La mitad del ancho del IC es:')
FCI    = nom_par.*FRCI


end