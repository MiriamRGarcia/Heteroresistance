%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate confidence intervals using Fisher Information Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FCI = Fisher_CI(sens_xT_FI, Cov, inv_Cov, sensCov, cl, Nexp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem sizes:
ntexp = size(sens_xT_FI, 1);
np    = size(sens_xT_FI, 2);

% ----------------------------------------------------------------------- %
% Calculate the Fisher Information Matrix:
FIM   = zeros(np, np);

for ir = 1:np
    sens_xT_pr = sens_xT_FI(:, ir);
    for ic = 1:np
        sens_xT_pc = sens_xT_FI(:, ic);
        FIM(ir, ic) = (sens_xT_pr.'*Cov*sens_xT_pc + 0.5*trace(inv_Cov*sensCov(:, :, ir)*inv_Cov*sensCov(:, :, ic)));
    end
end

% ----------------------------------------------------------%
% Inverse of Fisher information:

eigFIM   = eig(FIM)
detFIM   = prod(eigFIM)
rgFIM    = rank(FIM)

%pvar     = 1./sqrt(eigFIM);
%rank(FIM)
%det(FIM)
IFI  = inv(FIM)
pvar = abs(diag(IFI));

% ----------------------------------------------------------%
% Calculate the tStudent value at the required confidence level:
alph = 1 - cl;
pLo  = alph/2;
pUp  = 1 - alph/2;
ndf  = ntexp*Nexp - np;
tStud = tinv([pLo pUp], ndf);

% ----------------------------------------------------------%
% Calculate confidence interval:
FCI = tStud(2)*sqrt(pvar)

end