%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate confidence intervals using Fisher Information Matrix
% for normal iid error with constant covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FCI = Fisher_CI(Nexp, ntexp, par, nom_par, SensMatrix, sd, cl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT:
%%% SS = Sensitivity matrix,
%%% CV = Matrix with measurement covariances (constant),
%%% cl = Confidence level for the estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np = numel(par);

% ----------------------------------------------------------%
% Scale the sensitivity matrix by the measurement covariances:
inv_SD = diag(1/sd*ones(ntexp*Nexp, 1));
SS     = inv_SD*SensMatrix;

% ----------------------------------------------------------%
% Singular value decomposition of SS:
[~, SV, VV] = svd(SS);

dSV  = diag(SV)
IdSV = 1./(dSV.^2)
ISV  = diag(IdSV);

% ----------------------------------------------------------%
% Fisher Information and inverse:
FIM  = VV*(SV.'*SV)*VV.';
IFIM = VV*ISV*VV.';
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
nmFIM = cond(FIM)

fprintf('\n >> El rango de la Fisher es:')
rgFIM = rank(SV)

fprintf('\n >> El determinante de la Fisher es:')
detFIM = prod(dSV)

fprintf('\n >> La mitad del ancho del IC (OJO, RELATIVO!!) es:')
FRCI = tStud(2)*sqrt(pvar)

fprintf('\n >> La mitad del ancho del IC es:')
FCI  = nom_par.*FRCI

end