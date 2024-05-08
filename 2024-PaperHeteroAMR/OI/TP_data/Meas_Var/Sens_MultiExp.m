%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate average, variance and their sensitivities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xT, varxT, sens_xT, sens_varxT] = Sens_MultiExp(tmod, rr, Cexp, par, x0, ODEoptions)

% ----------------------------------------------------------------------- %
% Problem sizes:
Nexp = numel(Cexp);                                                        % Number of experiments,
nt   = numel(tmod);                                                        % Number of model times,
nr   = numel(rr);                                                          % Number of AMR degrees,
np   = numel(par);                                                         % Number of parameters,

% ----------------------------------------------------------------------- %
% Obtain parameter values:
mugS    = par(1);
mugR    = par(2);
alphg   = par(3);

mukmaxS = par(4);
bet     = par(5);
alphk   = par(6);

EC50k   = par(7);
Hk      = par(8);

xiSR    = par(9);
kxi     = par(10);

X0      = par(11);
lamb_IC = par(12);

% ----------------------------------------------------------------------- %
% Initial condition of the sensitivity system (same for each experiment):

% Initial condition for average and sensitivities:
s0 = [x0;
      zeros(nr*(np - 2),1)];                                               % Initial conditions are independent of parameters
                                                                           % unless for X0 and lamb_IC
% Ampliate with the initial condition for sensitivity to IC parameters:
f0_aux     = exp(-lamb_IC*rr);
s0_X0      = f0_aux/sum(f0_aux);
RR_aux     = repmat(rr.', nr, 1) - repmat(rr, 1, nr);    
s0_lamb_IC = X0/(sum(f0_aux)^2)*f0_aux.*(RR_aux*f0_aux);


s0 = [s0;
      s0_X0;
      s0_lamb_IC];

% Ampliate with the initial condition for variances and sensitivities:
s0 = [s0;
      zeros(1 + np, 1)];

% ----------------------------------------------------------------------- %
% Common calculations for ODEs (independent of the experiment):

% Calculate mutation-rates matrix:
RR  = repmat(rr, 1, nr) - repmat(rr.', nr, 1);                              % Auxiliary matrix to solve ODEs,
RR  = RR - triu(RR) + tril(RR).';

Xi = xiSR*exp(kxi*(1 - RR));
Xi = Xi - diag(diag(Xi));

% Auxiliary coefficient matrix:
AA_aux = Xi - diag(sum(Xi, 2));

% Derivative of the dynamics with respect to mugS:
d_mugS(1)  = 1;
d_mugS(nr) = 0;
d_mugS(2:(nr-1)) = mugR^2*(1 - rr(2:(nr-1)).^alphg)./(rr(2:(nr-1)).^alphg*(mugR - mugS) - mugR).^2;

AA_mugS    = diag(d_mugS);

% Derivative with respect to mugR:
d_mugR(1)  = 0;
d_mugR(nr) = 1;
d_mugR(2:(nr-1)) = mugS^2*rr(2:(nr-1)).^alphg./(mugR + rr(2:(nr-1)).^alphg*(mugS - mugR)).^2;

AA_mugR    = diag(d_mugR);

% Derivative with respect to alphg:
d_alphg(1)  = 0;
d_alphg(nr) = 0;
d_alphg(2:(nr-1)) = mugR*mugS*(mugR - mugS)*log(rr(2:(nr-1))).*rr(2:(nr-1)).^alphg./(mugR + rr(2:(nr-1)).^alphg*(mugS - mugR)).^2;

AA_alphg    = diag(d_alphg);

% Derivative with respect to xiSR:
d_xiSR  = exp(kxi*(1 - RR));
d_xiSR  = d_xiSR - diag(diag(d_xiSR));
AA_xiSR = d_xiSR - diag(sum(Xi/xiSR, 2));

% Derivative with respect to kxi:
d_kxi  = (1 - RR).*Xi;
d_kxi  = d_kxi - diag(diag(d_kxi));
AA_kxi = d_kxi - diag(sum(d_kxi, 2));

% Calculate growth rate:
mug = mugS*mugR./(mugR + (mugS - mugR)*rr.^alphg);

% Derivative of the growth rates with respect the different parameters:
% (constant in C):
der_mug = [d_mugS.' d_mugR.' d_alphg.' zeros(nr, np - 3)];

% Solve ODEs sensitivty system for each experiment:
xT         = zeros(nt, Nexp);
varxT      = zeros(nt, Nexp);
sens_xT    = zeros(nt, np, Nexp);
sens_varxT = zeros(nt, np, Nexp);

d_alphk    = zeros(nr,1);
d_bet      = zeros(nr,1);
d_Hk       = zeros(nr,1);
d_mukmaxS  = zeros(nr,1);
d_EC50k    = zeros(nr,1);

for iexp = 1:Nexp
    
    %-------------------------------------------------%
    % Calculate coefficient matrix of the state system:
    
    % Drug concentration:
    CC = Cexp(iexp);
    
    % Calculate kill rate at current time:
    mukmax = mukmaxS*bet^alphk*(1 - rr.^alphk)./(bet^alphk + rr.^alphk);
    HC     = CC^Hk/(CC^Hk + EC50k^Hk);
    muk    = mukmax*HC;

    AA     = AA_aux + diag(mug - muk);
    
    %-------------------------------------------------%
    % Calculate derivatives of the coefficient matrix:
    
     % Derivative of A with respect to mukmaxS:
    d_mukmaxS(1)  = 1;
    d_mukmaxS(nr) = 0;
    d_mukmaxS(2:(nr - 1)) = bet^alphk*(1 - rr(2:(nr - 1)).^alphk)./(bet^alphk + rr(2:(nr - 1)).^alphk);

    d_mukmaxS  = - HC*d_mukmaxS;

    AA_mukmaxS = diag(d_mukmaxS);

    % Derivative with respect to r50k:
    d_bet(1)  = 0;
    d_bet(nr) = 0;
    d_bet(2:(nr -1)) =  alphk*bet^(alphk - 1)*mukmaxS*(rr(2:(nr -1)).^alphk - 1).*rr(2:(nr -1)).^alphk./(bet^alphk + rr(2:(nr -1)).^alphk).^2;

    d_bet  = HC*d_bet;
    
    AA_r50k = diag(d_bet);

    % Derivative with respect to alph_k:
    d_alphk(1)  = 0;
    d_alphk(nr) = 0;
    d_alphk(2:(nr -1)) = bet^alphk*mukmaxS*rr(2:(nr -1)).^alphk.*((1 + bet^alphk)*log(rr(2:(nr - 1))) + (rr(2:(nr -1)).^alphk - 1)*log(bet))./(rr(2:(nr -1)).^alphk + bet^alphk).^2;

    d_alphk  = HC*d_alphk;

    AA_alphk = diag(d_alphk);

    % Derivative with respect to EC50k:
    d_EC50k(nr) = 0;
    d_EC50k(1:(nr - 1)) = - EC50k^(Hk - 1)*Hk*bet^alphk*mukmaxS*(rr(1:(nr -1)).^alphk - 1)./(rr(1:(nr -1)).^alphk + bet^alphk);
    d_EC50k = HC^2*d_EC50k;

    AA_EC50k = diag(d_EC50k);

    % Derivative with respect to Hk:
    if CC > 1e-20
        d_Hk(1:(nr -1)) = EC50k^Hk*bet^alphk*mukmaxS*(log(CC) - log(EC50k))*(rr(1:(nr -1)).^alphk - 1)./(rr(1:(nr -1)).^alphk + bet^alphk);
        d_Hk(nr) = 0;
        d_Hk = HC^2*d_Hk;
    else
        d_Hk(1:nr) = 0;
    end
    
    AA_Hk = diag(d_Hk);
    
    % Derivatives of kill rate:
    der_muk = - [zeros(nr, 3) d_mukmaxS d_bet d_alphk d_EC50k d_Hk zeros(nr, 4)];

    
    %-------------------------------------------------%
    % Calculate output sensitivities:
    [~, St_sens] = ode15s(@(t,s) Odes_Sens(t, s, AA, AA_mugS, AA_mugR, AA_alphg,...
                   AA_mukmaxS, AA_r50k, AA_alphk, AA_EC50k, AA_Hk, AA_xiSR, AA_kxi,...
                   mug, muk, der_mug, der_muk), tmod, s0, ODEoptions);
   
    % Total population:
    xx               = St_sens(1:nt, 1:nr);
    xT(1:nt, iexp)   = sum(xx, 2);
    
    % Sensitivities of the average:
    sensx_mugS    = St_sens(1:nt,nr + 1:2*nr);
    sensx_mugR    = St_sens(1:nt,2*nr + 1:3*nr);
    sensx_alphg   = St_sens(1:nt,3*nr + 1:4*nr);
    sensx_mukmaxS = St_sens(1:nt,4*nr + 1:5*nr);
    sensx_bet     = St_sens(1:nt,5*nr + 1:6*nr);
    sensx_alphk   = St_sens(1:nt,6*nr + 1:7*nr);
    sensx_EC50k   = St_sens(1:nt,7*nr + 1:8*nr);
    sensx_Hk      = St_sens(1:nt,8*nr + 1:9*nr);
    sensx_xiSR    = St_sens(1:nt,9*nr + 1:10*nr);
    sensx_kxi     = St_sens(1:nt,10*nr + 1:11*nr);
    sensx_X0      = St_sens(1:nt,11*nr + 1:12*nr);
    sensx_lambIC  = St_sens(1:nt,12*nr + 1:13*nr);
    
    sens_xT(1:nt, 1, iexp)  = sum(sensx_mugS, 2);
    sens_xT(1:nt, 2, iexp)  = sum(sensx_mugR, 2);
    sens_xT(1:nt, 3, iexp)  = sum(sensx_alphg, 2);
    sens_xT(1:nt, 4, iexp)  = sum(sensx_mukmaxS, 2);
    sens_xT(1:nt, 5, iexp)  = sum(sensx_bet, 2);
    sens_xT(1:nt, 6, iexp)  = sum(sensx_alphk, 2);
    sens_xT(1:nt, 7, iexp)  = sum(sensx_EC50k, 2);
    sens_xT(1:nt, 8, iexp)  = sum(sensx_Hk, 2);
    sens_xT(1:nt, 9, iexp)  = sum(sensx_xiSR, 2);
    sens_xT(1:nt, 10, iexp) = sum(sensx_kxi, 2);
    sens_xT(1:nt, 11, iexp) = sum(sensx_X0, 2);
    sens_xT(1:nt, 12, iexp) = sum(sensx_lambIC, 2);

    % Variance of the total population:
    varxT(1:nt, iexp) = St_sens(1:nt, 13*nr + 1);
    
    % Sensitivities of the variance of the total population:
    for ip = 1:np
        sens_varxT(1:nt, ip, iexp) = St_sens(1:nt, 13*nr + 1 + ip);
    end
    
end

          
end