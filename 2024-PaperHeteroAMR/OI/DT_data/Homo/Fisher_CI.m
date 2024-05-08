%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate confidence intervals using Fisher Information Matrix
% for normal iid error with constant covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FCI = Fisher_CI(Nexp, ntexp, par, nom_par, SensMatrix, sd, cl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------%
% Indices of the parameters to calculate FIM:
pind   = [1 2 3 4 5 6 7 8 9 10 11 12];

% ----------------------------------------------------------%
% Scale the sensitivity matrix by the measurement covariances:
inv_SD = diag((1/sd)*ones(ntexp*Nexp, 1));
SS     = inv_SD*SensMatrix(:, pind);
np     = numel(par(pind));

% ----------------------------------------------------------%
% Singular value decomposition of SS:
[~, SV, VV] = svd(SS);

dSV   = diag(SV);
SV2   = diag(dSV.^2);

IdSV2 = 1./diag(SV2);
ISV2  = diag(IdSV2);

% ----------------------------------------------------------%
% Calculate Fisher Information and inverse:
FIM  = VV*SV2*VV.';
IFIM = VV*ISV2*VV.';

% Comprobate if FIM is ill-conditioned:
rcondFIM = rcond(FIM);

fprintf('\n >> El rango de la Fisher es: %d', rank(SV))
fprintf('\n >> El determinante de la Fisher es: %d', prod(dSV))
fprintf('\n >> El numero de condicion de la Fisher es: %d', cond(FIM))
fprintf('\n >> El numero de condicion reciproco de la Fisher es: %d', rcondFIM)

% Regularise FIM if ill-conditioned:
if (rcondFIM > 1e-50) && (rcondFIM < 1e-10)
    cc    = 50;
    Rdelt = max(0, (dSV(1) - cc*dSV(np))/(cc - 1)); 
    FIM   = (FIM + Rdelt*eye(np))/(1 + Rdelt);
    
    rcondFIM = rcond(FIM);
    
    fprintf('\n >> El numero de condicion reciproco tras la regularizacion es: %d', rcondFIM)
    
    [auxUU, auxSV, auxVV] = svd(FIM);
    dauxSV                = diag(auxSV);
    IFIM                  = auxVV*diag(1./dauxSV)*auxUU.';
end

% Tykhonov regularisation (does not work):
% if (condIFIM > 1e-50) &&  (condIFIM < 1e-10)
%     alphTik = min(diag(FIM))*0.01; % regularization of the diagonals
%     FIM     = FIM + alphTik*((par(pind) - nom_par(pind)).'*eye(np)*(par(pind) - nom_par(pind)));
%     IFIM    = pinv(FIM);
% end

% Estimate of parameter variances:
pvar = diag(IFIM);

% ----------------------------------------------------------%
% Calculate the tStudent value at the required confidence level:
alph  = 1 - cl;
pLo   = alph/2;
pUp   = 1 - alph/2;
ndf   = ntexp*Nexp - np;
tStud = tinv([pLo pUp], ndf);

% ----------------------------------------------------------%
% Calculate confidence interval:
FRCI   = tStud(2)*sqrt(pvar);
fprintf('\n >> La mitad del ancho del IC (OJO, RELATIVO!!) es: [');
fprintf('%g ', FRCI);
fprintf(']');

FCI    = nom_par(pind).*FRCI;
fprintf('\n >> La mitad del ancho del IC es: [');
fprintf('%g ', FCI);
fprintf(']');



end