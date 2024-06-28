%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIM_ConfInt: Calculate FIM confidence intervals from sensitivities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FIM, FIMConfInt] = FIM_ConfInt(FIM_pars, pars_nom, pars_var,...
    SensMatrix, normSensMatrix, CovMatrix, confLev, noise, Nexp, ntexp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% FIM_pars       =
% pars_nom       =
% SensMatrix     =
% normSensMatrix =
% CovMatrix      = 
% confLev        =
% noise
% Nexp           =
% ntexp          =
%
% OUTPUT:
% FIM            = Fisher Information Matrix (scaled by nominalparameters);
% FIMConfInt     = FIM confidence intervals (size np x 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np = numel(FIM_pars);

% ----------------------------------------------------------------------- %
% Scaling the sensitivity matrix by the standard deviation matrix:
dgCov         = diag(CovMatrix);
dgSD          = sqrt(dgCov);
inv_SD        = diag(1./dgSD);
sc_SensMatrix = inv_SD*SensMatrix;

if strcmp(noise, 'MNHe')
    % Obtain parameters of variance:
    var_a = pars_var(1);
    var_b = pars_var(2);

    ndata = size(CovMatrix, 1);

    % ------------------------------------------------------------------- %
    % Calculate the term of FIM in the sensitivity matrix:

    % Ampliate with derivatives of the average to var_b:
    sc_SensMatrix = [sc_SensMatrix zeros(ndata, 1)];

    % Singular value decomposition:
    [~, SensMatrixSV, SensMatrix_VV] = svd(sc_SensMatrix);

    SensMatrixSV = diag(SensMatrixSV);     
    SV2          = diag(SensMatrixSV.^2); 

    % First term of FIM:
    FIM = SensMatrix_VV*SV2*SensMatrix_VV.';

    % ------------------------------------------------------------------- %
    % Calculate the term of FIM in covariances:
    % Note: this is true only taking the rescaled parameters p^* = p/ptrue
    aux   = (dgCov/var_a).^(1/var_b);
    sc_normSensMatrix  = [var_b*normSensMatrix var_b*log(aux)];

    % Singular value decomposition of the scaled sensitivity matrix:
    [~, normSensMatrixSV, normSensMatrix_VV] = svd(sc_normSensMatrix);

    dgnormSensMatrixSV = diag(normSensMatrixSV);    
    normSensMatrixSV2  = diag(dgnormSensMatrixSV.^2); 
    
    % Trace of covariance matrix:
    TrCov = normSensMatrix_VV*normSensMatrixSV2*normSensMatrix_VV.';

    % Term of the FIM in the traces of covariance matrix sensitivities:
    FIM = FIM + 0.5*TrCov;

    % ------------------------------------------------------------------- %
    % Calculate the inverse of FIM:

    % Cholesky factorisation of FIM:
    Chol_FIM = chol(FIM);

    % Inverse of the FIM from svd of RR_FIM:
    [~, Chol_FIM_SV, Chol_FIM_V] = svd(Chol_FIM);

    dChol_FIM_SV  = diag(Chol_FIM_SV);         
    Chol_FIM_ISV2 = diag(1./(dChol_FIM_SV.^2)); 

    IFIM = Chol_FIM_V*Chol_FIM_ISV2*Chol_FIM_V.';  
    
else
    
    % ------------------------------------------------------------------- %
    % Singular value decomposition of Sensitivity Matrix:
    [~, SensMatrixSV, SensMatrix_VV] = svd(sc_SensMatrix);

    SensMatrixSV = diag(SensMatrixSV);
    dgSV2        = SensMatrixSV.^2;
    SV2          = diag(dgSV2);
    InvdgSV2     = 1./dgSV2;
    InvSV2       = diag(InvdgSV2);

    % ------------------------------------------------------------------- %
    % Calculate FIM and inverse FIM:
    FIM  = SensMatrix_VV*SV2*SensMatrix_VV.';
    IFIM = SensMatrix_VV*InvSV2*SensMatrix_VV.';

end

% ----------------------------------------------------------%
% Calculate the tStudent value at the required confidence level:
alph  = 1 - confLev;
pLo   = alph/2;
pUp   = 1 - alph/2;
ndf   = ntexp*Nexp - np;
tStud = tinv([pLo pUp], ndf);

% ----------------------------------------------------------%
% Calculate confidence interval:
pvar = diag(IFIM);

% Scale by nominal parameters:
RFIMConfInt = tStud(2)*sqrt(pvar);
FIMConfInt  = [FIM_pars-pars_nom.*RFIMConfInt FIM_pars+pars_nom.*RFIMConfInt];



end