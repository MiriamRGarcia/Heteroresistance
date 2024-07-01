%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIM_CI: Calculate FIM confidence intervals from sensitivities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FIM, CI] = FIM_CI(pars_opt, pars_var, SensMatrix, normSensMatrix,...
                            CovMatrix, confLev, noise, m_r, m_texp, m_e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% OUTPUT:
% FIM = Fisher Information Matrix (with log-scaling);
% CI  = FIM-based confidence intervals (size: m_p x 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem sizes:
m_p = numel(pars_opt);

% ----------------------------------------------------------------------- %
% Scaling the sensitivity matrix by the standard deviations:
dgCov         = diag(CovMatrix);
dgSD          = sqrt(dgCov);
inv_SD        = diag(1./dgSD);
sc_SensMatrix = inv_SD*SensMatrix;

if strcmp(noise, 'MNHe')
    
    % Consider var_b in the number of parameters:
    m_p_aux = m_p + 1;
    
    % Obtain parameters for variance:
    var_a_opt = pars_var(1);
    var_b_opt = pars_var(2);
    
    % Size of the dataset:
    m_d = size(CovMatrix, 1);

    % ------------------------------------------------------------------- %
    % Calculate the term of FIM in the sensitivity matrix:

    % Ampliate with derivatives of the average to var_b:
    sc_SensMatrix = [sc_SensMatrix zeros(m_d, 1)];

    % Singular value decomposition:
    [~, SensMatrixSV, SensMatrix_VV] = svd(sc_SensMatrix);

    SensMatrixSV = diag(SensMatrixSV);     
    SV2          = diag(SensMatrixSV.^2); 

    % First term of FIM:
    FIM = SensMatrix_VV*SV2*SensMatrix_VV.';

    % ------------------------------------------------------------------- %
    % Calculate the term of FIM in covariances:
    % Note: this is true only taking the rescaled parameters:
    % pars^* = log(pars):
    aux   = (dgCov/var_a_opt).^(1/var_b_opt);
    sc_normSensMatrix  = [var_b_opt*normSensMatrix var_b_opt*log(aux)];

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
    
    pars_opt = [pars_opt;pars_var(end)];
    
else
    
    m_p_aux = m_p;
    
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
ndf   = m_texp*m_e - m_p_aux;
tStud = tinv([pLo pUp], ndf);

% ----------------------------------------------------------%
% Calculate FIM-based confidence interval:
pvar = diag(IFIM);

% Half amplitude of the confidence interval 
% for reescaled parameters, pars^* = log(pars):
RFIMConfInt = tStud(2)*sqrt(pvar);

% Half amplitude of the confidence interval 
% for model parameters, pars:
FCI = pars_opt.*RFIMConfInt;

% Confidence interval:
CI  = [pars_opt - FCI, pars_opt + FCI];

% ----------------------------------------------------------------------- %
% Print results:
fprintf('\n>> The determinant of FIM is: %.4e', det(FIM))
fprintf('\n>> The reciprocal condition number of FIM is: %.4e', rcond(FIM))
fprintf('\n>> The FIM confidence intervals for the model parameters are:');

if m_r > 2
    par_names = {'b_S','b_R','alpha_b','d_Smax','alpha_d','beta_d',...
                 'EC_50d','H_d','xi_SR', 'k_xi','N_T0','lambda_T0'};
else
    par_names = {'b_S','b_R','d_Smax','EC_50d','H_d','xi_SR','N_T0','lambda_T0'};
end

for ip = 1:m_p
    aux_par_names = cell2mat(par_names(:,ip));
    fprintf('\n>> %s ---> (%.4e - %.4e, %.4e + %.4e)', aux_par_names, ...
                             pars_opt(ip), FCI(ip), pars_opt(ip), FCI(ip))
end

if strcmp(noise, 'MNHe')
    fprintf('\n>> var_b ---> (%.4e - %.4e, %.4e + %.4e)\n', pars_opt(m_p_aux), FCI(m_p_aux), pars_opt(m_p_aux), FCI(m_p_aux))
end

end