%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate confidence intervals using Fisher Information Matrix for iid 
% normal measurement noise and known variance + deterministic data
% using CFUS/mL as measured variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

% ----------------------------------------------------------------------- %
% Experimental setting:

% Parameter values:
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

Y0      = 6;
lamb_IC = 50;

par     = [mugS;mugR;alphg;mukmaxS;bet;alphk;EC50k;Hk;xiSR;kxi;Y0;lamb_IC];
nom_par = par;
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

% Initial condition:
f0  = exp(-lamb_IC*rr);

y0  = log(10)*Y0 - lamb_IC*rr - log(sum(f0))*ones(nr, 1);

% ----------------------------------------------------------------------- %
% Calculate average and sensitivities:
ODEoptions = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);

[yT, sens_yT] = Sens_MultiExp(tmod, rr, Cexp, par, y0, ODEoptions);

% Log-10 measures of the average:
yT = log10(xT);

% Sensitivities of the log10:
sens_yT = zeros(nt, np, Nexp);

for ip = 1:np
    aux = reshape(sens_xT(1:nt, ip, 1:Nexp), nt, Nexp);
    aux = aux./xT/log(10);
    sens_yT(1:nt, ip, 1:Nexp) = aux;
end

% ----------------------------------------------------------------------- %
% (Optional) Plot sensitivities and normalised sensitivities if desired:
figs = 0;

if figs
    Plot_Sens(tmod, Cexp, par, yT, sens_yT)
end


% ----------------------------------------------------------------------- %
% Calculate sensitivity matrix:
SensMatrix = [];
for ip = 1:np
    aux = reshape(sens_yT(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    
    % Scaling by the parameter:
    %aux = nom_par(ip)*aux;
    aux = par(ip)*log(10)*aux;
    
    SensMatrix = [SensMatrix reshape(aux.',[],1)];
end

% Multi-exp normalised sensitivity matrix using the sampling times:       
scSensMatrix = [];
for ip = 1:np
    aux = reshape(sens_yT(texp_ind, ip, 1:Nexp), ntexp, Nexp);
    aux = par(ip)*aux./yT(texp_ind, 1:Nexp);
    scSensMatrix = [scSensMatrix reshape(aux.',[],1)];
end

[~,SVS, VS] = svd(SensMatrix);
[~,SVscS, VscS] = svd(scSensMatrix);

rgS   = rank(SVS);
rgscS = rank(SVscS);

% ----------------------------------------------------------------------- %
% Calculate correlation matrix:
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
% Calculate confidence intervals from Fisher Information Matrix:
cl  = 0.95;                                                                % Confidence level,
sd  = 0.5;                                                                 % Constant standard deviation of measures,

FCI = Fisher_CI(Nexp, ntexp, par, nom_par, SensMatrix, sd, cl);                     % Half-lenght of the confidence interval for each parameter,